####!/usr/bin/env python2.7

###
###   This script is part of Sequenza
###   http://www.cbs.dtu.dk/biotools/sequenza/
###

import argparse, os, sys, gzip, math, multiprocessing, time, logging, json, gc
from itertools import izip_longest
from multiprocessing.pool import ThreadPool
from functools import partial
from multiprocessing.queues import SimpleQueue
from contextlib import contextmanager
from shutil     import rmtree
from subprocess import Popen, check_call, PIPE
from tempfile   import mkdtemp

VERSION = "2.1.0"
DATE    = "07 October 2014"
AUTHOR  = "Favero Francesco"
MAIL    = "favero@cbs.dtu.dk"

@contextmanager
def named_pipe():
   '''
   Thanks to google for this nice solution :)
   '''
   dirname = mkdtemp()
   try:
      path = os.path.join(dirname, "temp.fifo")
      os.mkfifo(path)
      yield path
   finally:
      rmtree(dirname)

def check_dir(directory):
   '''
   Check the given directory and return the whole path
   or set the default to working directory.
   '''
   if directory == '':
      directory = None
   if directory == None:
      directory = os.getcwd()
   if os.path.isdir(directory):
      directory = os.path.abspath(directory)
   else:
      sys.exit("Not existing directory, please create the directory first")
   return directory

def walklevel(directory, level=1):
   '''

   '''
   directory = directory.rstrip(os.path.sep)
   assert os.path.isdir(directory)
   num_sep = directory.count(os.path.sep)
   for root, dirs, files in os.walk(directory):
      yield root, dirs, files
      num_sep_this = root.count(os.path.sep)
      if num_sep + level <= num_sep_this:
         del dirs[:]

def xopen(filename, mode='r'):
   '''
   Replacement for the "open" function that can also open
   files that have been compressed with gzip. If the filename ends with .gz,
   the file is opened with gzip.open(). If it doesn't, the regular open()
   is used. If the filename is '-', standard output (mode 'w') or input
   (mode 'r') is returned.
   '''
   assert isinstance(filename, basestring)
   if filename == '-':
      return sys.stdin if 'rb' in mode else sys.stdout
   if filename.endswith('.gz'):
      return gzip.open(filename, mode)
   else:
      return open(filename, mode)

class IterableQueue(SimpleQueue):
   _sentinel = object()
   def next(self):
      item = self.get()
      if item == "EOF":
         self.close()
         raise StopIteration
      return item
   def __iter__(self):
      return (iter(self.next, self._sentinel))
   def close(self):
      self.put(self._sentinel)


def grouper(n, iterable, padvalue=None):
   """grouper(3, 'abcdefg', 'x') -->
   ('a','b','c'), ('d','e','f'), ('g','x','x')"""
   return izip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

def pileup_partial_split(pileup_line):
   '''
   split pileup line in chromosome, position and rest of line,
   without even strip
   '''
   chromosome, position, line = pileup_line.split('\t', 2)
   return [chromosome, int(position), line]

class multiPileups:
   '''
   Define an iterable object merging two pileups
   '''
   def __init__(self, p1, p2):
      self.p1 = p1
      self.p2 = p2
      self._last_chromosome = None
   _sentinel        = object()
   def next(self):
      try:
         chromosome1, position1, line1 = pileup_partial_split(self.p1.next())
         chromosome2, position2, line2 = pileup_partial_split(self.p2.next())
      except ValueError:
         raise StopIteration
      if chromosome1 == chromosome2:
         self._last_chromosome = chromosome1
      going_on = True
      while going_on:
         '''
         If one of the files finishes, it will raise StopIteration
         and will stop the iteration cleanly
         '''
         try:
            if chromosome1 == chromosome2:
               self._last_chromosome = chromosome1
               if position1 == position2:
                  ref1, line1 = line1.split('\t', 1)
                  ref2, line2 = line2.split('\t', 1)
                  return [chromosome1, position1, ref1, line1.strip(), line2.strip()]
                  going_on = False
               elif position1 > position2:
                  chromosome2, position2, line2 = pileup_partial_split(self.p2.next())
               elif position1 < position2:
                  chromosome1, position1, line1 = pileup_partial_split(self.p1.next())
            else:
               if chromosome1 != self._last_chromosome:
                  chromosome2, position2, line2 = pileup_partial_split(self.p2.next())
               else:
                  chromosome1, position1, line1 = pileup_partial_split(self.p1.next())
         except StopIteration:
            going_on = False
            raise StopIteration
   def close(self):
      self.p1.close()
      self.p2.close()
   def __iter__(self):
      return (iter(self.next, self._sentinel))

def next_gcfile(self):
   gc_line = self.gc.next()
   try:
      pos, self._last_GC    = gc_line.split("\t", 1)
      self._last_window_s   = int(pos)
      self._last_window_e   = self._last_window_s + self._last_span
      return self
   except ValueError:
      line, chrom, span     = gc_line.strip().split()
      self._chromosome = chrom.split('chrom=')[1]
      self._last_span       = int(span.split('span=')[1])
      pos, self._last_GC    = self.gc.next().strip().split("\t", 1)
      self._last_window_s   = int(pos)
      self._last_window_e   = self._last_window_s + self._last_span
      return self
   except ValueError:
      if self._chromosome:
         raise StopIteration
      else:
         sys.exit("Error: the GC-content file is not in the expected format. See for example http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/gc5Base/ ")

class GCmultiPileups:
   '''
   Add GC-content to the pileup, and still make an iterable object
   '''
   def __init__(self,multiPileups,gc):
      self.gc          = gc
      self             = next_gcfile(self)
      self.mpileup     = multiPileups
      self._last_chromosome = None
   _sentinel        = object()
   def next(self):
      self.pile_line   = self.mpileup.next()
      going_on = True
      while going_on:
         if self._chromosome == self.pile_line[0]:
            self._last_chromosome = self._chromosome
            if self.pile_line[1] >= self._last_window_s and self.pile_line[1] < self._last_window_e:
               self.pile_line.append(self._last_GC.strip())
               return self.pile_line
               going_on = False
            elif self.pile_line[1] < self._last_window_s:
               self.pile_line = self.mpileup.next()
            elif self.pile_line[1] >= self._last_window_e:
               self = next_gcfile(self)
         else:
            if self._last_chromosome != self._chromosome and self._last_chromosome != None:
               self.pile_line = self.mpileup.next()
            else:
               self = next_gcfile(self)
   def close(self):
      self.mpileup.close()
      self.gc.close()
   def __iter__(self):
      return (iter(self.next, self._sentinel))

class GCmultiPileupsAltDepth:
   '''
   Add alternative depth to the GCpileup, and again, keep it iterable
   '''
   def __init__(self,GCmultiPileups,p2):
      self.mpileup = GCmultiPileups
      self.p2      = p2
      self._last_chromosome = None
   _sentinel        = object()
   def next(self):
      self.pile_line = self.mpileup.next()
      try:
         chromosome2, position2, line2 = pileup_partial_split(self.p2.next())
      except ValueError:
         raise StopIteration
      if self.pile_line[0] == chromosome2:
         self._last_chromosome = self.pile_line[0]
      going_on = True
      while going_on:
         '''
         If one of the files ever finish it would rise StopIteration
         and it would stop the iteration cleanly
         '''
         try:
            if self.pile_line[0] == chromosome2:
               self._last_chromosome = self.pile_line[0]
               if self.pile_line[1] == position2:
                  ref2, depth2, line2 = line2.split('\t', 2)
                  self.pile_line.append(depth2.strip())
                  return self.pile_line
                  going_on = False
               elif self.pile_line[1] > position2:
                  chromosome2, position2, line2 = pileup_partial_split(self.p2.next())
               elif self.pile_line[1] < position2:
                  self.pile_line = self.mpileup.next()
            else:
               if self.pile_line[0] != self._last_chromosome:
                  chromosome2, position2, line2 = pileup_partial_split(self.p2.next())
               else:
                  self.pile_line = self.mpileup.next()
         except StopIteration:
            going_on = False
            raise StopIteration
   def close(self):
      self.mpileup.close()
      self.p2.close()
   def __iter__(self):
      return (iter(self.next, self._sentinel))

def parse_pileup(ref_base, line, qlimit=20, qformat='sanger'):
   '''
   Parse the pileup format
   '''
   depth, mut_list, mut_qual = line
   ref_base = ref_base.upper()
   depth    = int(depth)
   if ref_base != "N":
      freq_dict = parse_pileup_seq(mut_list, mut_qual, depth, ref_base, qlimit, qformat)
      return [ref_base, depth, freq_dict['A'], freq_dict['C'], freq_dict['G'], freq_dict['T'], freq_dict['Z']]
   else:
      return [ref_base, depth, 0, 0, 0, 0, [0, 0, 0, 0]]

def parse_pileup_str(line, min_depth, qlimit=20, qformat='sanger'):
   '''
   Parse the pileup format
   '''
   if line:
      line   = line.strip()
      try:
         chromosome, n_base, ref_base, depth, mut_list, mut_qual = line.split()
         rd     = int(depth)
         rebase = ref_base.upper()
         if rd >= min_depth and rebase != "N":
            freq_dict = parse_pileup_seq(mut_list, mut_qual, rd, rebase, qlimit, qformat)
            line = [chromosome, n_base, rebase, depth, freq_dict['A'], freq_dict['C'], freq_dict['G'], freq_dict['T'], ':'.join(map(str, freq_dict['Z']))]
            return '\t'.join(map(str,line))
         else:
            pass
      except ValueError:
         pass
   else:
      pass

def parse_pileup_seq(seq, quality, depth, reference, qlimit=20, qformat='sanger', noend=False, nostart=False):
   '''
   Parse the piled up sequence and return the frequence of mutation.
   '''
   nucleot_dict = {'A':0, 'C':0, 'G':0, 'T':0}
   strand_dict  = {'A':0, 'C':0, 'G':0, 'T':0}
   reference    = reference.strip().upper()
   H            = 0
   n            = 0
   if qformat == 'sanger':
      qlimit = qlimit + 33
   elif qformat == 'illumina':
      qlimit = qlimit + 64
   else:
      sys.exit("Supported quality format are only \"illumina\" and \"sanger\"(default).")
   quality      = quality.strip()
   block = {'seq':'','length':0}
   start     = False
   del_ins   = False
   l_del_ins = ''
   last_base = ''
   ins_del_length = 0
   for base in seq.strip():
      if block['length'] == 0:
         if base == '$':
            if noend:
               nucleot_dict[last_base] -= 1
         elif base == '^':
            start            = True
            block['length'] += 1
            block['seq']     = base
         elif base == '+' or base == '-':
            del_ins          = True
            block['length'] += 1
            block['seq']     = base
         elif base == '.' or base == ',':
            if ord(quality[n]) >= qlimit:
               nucleot_dict[reference] += 1
               if base == '.': strand_dict[reference] += 1
            last_base        = reference
            n += 1
         elif base.strip().upper() in nucleot_dict:
            if ord(quality[n]) >= qlimit:
               nucleot_dict[base.strip().upper()] += 1
               if base.strip().isupper(): strand_dict[base.strip().upper()] += 1
            last_base        = base.strip().upper()
            n +=1
         else:
            n += 1
      else:
         if start:
            block['length'] += 1
            block['seq']    += base
            if block['length'] == 3:
               if not nostart:
                  if base == '.' or base == ',':
                     if ord(quality[n]) >= qlimit:
                        nucleot_dict[reference] += 1
                        if base == '.': strand_dict[reference] += 1
                  elif base.strip().upper() in nucleot_dict:
                     if ord(quality[n]) >= qlimit:
                        nucleot_dict[base.strip().upper()] += 1
                        if base.strip().isupper(): strand_dict[base.strip().upper()] += 1
               block['length'] = 0
               block['seq']    = ''
               start           = False
               n += 1
         elif del_ins:
            if base.isdigit():
               l_del_ins       += base
               block['seq']    += base
               block['length'] += 1
            else:
               ins_del_length = int(l_del_ins)+ 1 + len(l_del_ins)
               block['seq']    += base
               block['length'] += 1
               if block['length'] == ins_del_length:
                  block['length'] = 0
                  block['seq']    = ''
                  l_del_ins       = ''
                  ins_del         = False
                  ins_del_length  = 0
      # Debug line
      #print str(n) + " " + base + " " + seq
   if n == depth:
      nucleot_dict["Z"] = [strand_dict['A'], strand_dict['C'], strand_dict['G'], strand_dict['T']]
      return nucleot_dict
   else:
      #print seq + " " + str(depth)+" " + str(n)
      sys.exit('Something went wrong on the block counting!!')


def line_worker(line, depth_sum, qlimit=20, qformat='sanger', hom_t=0.85, het_t=0.35, alt_pileup=False):
   '''
   After the 3 files are syncronized we need to transform the pileup in a
   readable format, find the alleles, and compute the allele frequency in tumor
   '''
   if line:
      bases_list = ['A', 'C', 'G', 'T']
      if alt_pileup:
         chromosome, position, ref, p1_str, p2_str, gc, alt_depth = line
      else:
         chromosome, position, ref, p1_str, p2_str, gc = line
      p1_list = p1_str.split()
      p2_list = p2_str.split()
      if int(p1_list[0]) + int(p2_list[0]) >= depth_sum and len(p1_list) == len(p2_list):
         p1_mu  = parse_pileup(ref, p1_list, qlimit, qformat)
         sum_p1 = float(sum(p1_mu[2:6]))
         if sum_p1 > 0:
            p1_freq = [x/sum_p1 for x in p1_mu[2:6]]
            sort_freq_p1 = sorted(p1_freq, reverse = True)
            if sort_freq_p1[0] >= hom_t:
               # Homozygous positions here
               i        = p1_freq.index(sort_freq_p1[0])
               p2_mu    = parse_pileup(ref, p2_list, qlimit, qformat)
               sum_p2   = float(sum(p2_mu[2:6]))
               if sum_p2 > 0:
                  p2_freq  = [x/sum_p2 for x in p2_mu[2:6]]
                  homoz_p2 = p2_freq[i]
                  no_zero_idx   = [val for val in range(len(p2_freq)) if p2_freq[val] > 0]
                  no_zero_bases = [str(bases_list[ll])+str(round(p2_freq[ll], 3)) for ll in no_zero_idx if ll != i]
                  if no_zero_bases == []:
                     no_zero_bases = '.'
                     strands_bases = '0'
                  else:
                     no_zero_bases = ":".join(map(str, no_zero_bases))
                     strands_bases = [str(bases_list[ll])+str(round(p2_mu[6][ll]/float(p2_mu[2+ll]), 3)) for ll in no_zero_idx if ll != i]
                     strands_bases = ":".join(map(str, strands_bases))
                  #homoz_p2 = p2_mu[2 + i]/sum_p2
                  # chromosome, n_base, base_ref, depth_normal, depth_sample, depth.ratio, Af, Bf, zygosity.normal, GC-content, reads above quality, AB.ref, AB.tum, percentage AB.tum in fw
                  if alt_pileup:
                     line_out = [chromosome, position, p1_mu[0], alt_depth, p2_mu[1], round(p2_mu[1]/float(alt_depth), 3), round(homoz_p2, 3), 0, 'hom', gc, int(sum_p2), bases_list[i], no_zero_bases, strands_bases]
                  else:
                     line_out = [chromosome, position, p1_mu[0], p1_mu[1], p2_mu[1], round(p2_mu[1]/float(p1_mu[1]), 3), round(homoz_p2, 3), 0, 'hom', gc, int(sum_p2), bases_list[i], no_zero_bases, strands_bases]
                  return line_out
               else:
                  pass
            else:
               if sort_freq_p1[1] >= het_t:
                  # Heterozygous position here
                  allele = list()
                  for b in xrange(4):
                     if p1_freq[b] >= het_t:
                        allele.append(bases_list[b])
                  if len(allele) == 2:
                     p2_mu    = parse_pileup(ref, p2_list, qlimit, qformat)
                     sum_p2   = float(sum(p2_mu[2:6]))
                     if sum_p2 > 0:
                        i        = bases_list.index(allele[0])
                        ii       = bases_list.index(allele[1])
                        het_a_p2 = p2_mu[2 + i]/sum_p2
                        het_b_p2 = p2_mu[2 + ii]/sum_p2
                        if  het_a_p2 >= het_b_p2:
                           if alt_pileup:
                              line_out = [chromosome, position, p1_mu[0], alt_depth, p2_mu[1], round(p2_mu[1]/float(alt_depth), 3), round(het_a_p2 , 3), round(het_b_p2 , 3) , 'het', gc, int(sum_p2), bases_list[i]+bases_list[ii], ".", "0"]
                           else:
                              line_out = [chromosome, position, p1_mu[0], p1_mu[1], p2_mu[1], round(p2_mu[1]/float(p1_mu[1]), 3), round(het_a_p2 , 3), round(het_b_p2 , 3) , 'het', gc, int(sum_p2), bases_list[i]+bases_list[ii], ".", "0"]
                           return line_out
                        elif het_a_p2 < het_b_p2:
                           if alt_pileup:
                              line_out = [chromosome, position, p1_mu[0], alt_depth, p2_mu[1], round(p2_mu[1]/float(alt_depth), 3), round(het_b_p2 , 3), round(het_a_p2 , 3) , 'het', gc, int(sum_p2), bases_list[ii]+bases_list[i], ".", "0"]
                           else:
                              line_out = [chromosome, position, p1_mu[0], p1_mu[1], p2_mu[1], round(p2_mu[1]/float(p1_mu[1]), 3), round(het_b_p2 , 3), round(het_a_p2 , 3) , 'het', gc, int(sum_p2), bases_list[ii]+bases_list[i], ".", "0"]
                           return line_out
                  else:
                     pass
         else:
            pass
      else:
         pass
   else:
      pass

class abfreReduce:
   def __init__(self, abf, w):
      self.abf = abf
      self.w   = w
      self._status = 1
      self._header = next(self.abf).strip()
      line    = next(self.abf).strip()
      line_ls = line.split('\t')
      self.__refresh__(line_ls, line)
   _sentinel        = object()
   def next(self):
      if self._status == 1:
         try:
            line    = next(self.abf).strip()
            line_ls = line.split('\t')
         except StopIteration:
            self._status = 0
            self.__do_dict__()
            #return [line_ls[1], self.w + self._last_edge, self._n]
            return self.line_dict    
         if self.w + self._last_edge < int(line_ls[1]) or self._last_chromosome != line_ls[0]:
            self.__do_dict__()
            line_dict = self.line_dict
            self.__refresh__(line_ls, line)
            #return [line_ls[1], self.w + self._last_edge, self._n]
            return line_dict 
         else:
            self.__addline__(line_ls, line)
      elif self._status == 0:
         raise StopIteration
   def __refresh__(self, line_ls, line):
      self._last_chromosome = line_ls[0]
      self._last_edge       = int(line_ls[1])
      self._last_position   = self._last_edge
      self._n_depth         = int(line_ls[3])
      self._t_depth         = int(line_ls[4])
      self._ratio           = float(line_ls[5])
      self._gc              = float(line_ls[9])
      self._n               = 1
      self.line_dict = {'top': '', 'middle': [], 'end': ''}
      if line_ls[12]+line_ls[8] != '.hom':
         self.line_dict['top'] = line_ls
   def __addline__(self, line_ls, line):
      self._last_position    = int(line_ls[1])
      self._n_depth         += int(line_ls[3])
      self._t_depth         += int(line_ls[4])
      self._ratio           += float(line_ls[5])
      self._gc              += float(line_ls[9])
      self._n               += 1
      if line_ls[12]+line_ls[8] != '.hom':
         self.line_dict['middle'].append(line_ls)
   def __do_dict__(self):
      gc = str(int(round(self._gc/self._n, 0)))
      avg_line = [self._last_chromosome, self._last_position, 'N', self._n_depth/self._n,
                  self._t_depth/self._n, round(self._ratio/self._n, 3), 1.0, 0, 'hom',
                  gc , self._n, 'N', '.', '0']
      self.line_dict['end'] = map(str, avg_line)
      if self.line_dict['top'] == '':
         avg_line[1] = self._last_edge
         self.line_dict['top'] = map(str, avg_line)
      else:
         self.line_dict['top'][9] = gc
      for mid in self.line_dict['middle']:
         mid[9] = gc
   def close(self):
      self.abf.close()
   def __iter__(self):
      return (iter(self.next, self._sentinel))

def stream_fasta(fa_file, window_size, out_queue):
   """
   Read and stream a fasta file in a Pipe
   """
   tmp_line   = ''
   for line in fa_file:
      line = line.strip()
      if line.startswith('>'):
         chromosome = line
         if tmp_line == '':
            #print chromosome
            out_queue.put(chromosome)
         else:
            groups = grouper(window_size, tmp_line.upper())
            for group in groups:
               #print group
               out_queue.put(group)
            tmp_line = ''
            #print chromosome
            out_queue.put(chromosome)
      else:
         tmp_line = tmp_line + line
         if len(tmp_line)%window_size == 0:
            groups = grouper(window_size, tmp_line.upper())
            for group in groups:
               #print group
               out_queue.put(group)
            tmp_line = ''
         else:
            #print tmp_line
            pass
   if tmp_line == '':
      pass
   else:
      groups = grouper(window_size, tmp_line.upper())
      for group in groups:
         #print group
         out_queue.put(group)
   out_queue.put('EOF')

def seq_map(seq_list):
   """
   ...
   """
   counter = {"A":0, "C":0, "G":0, "N":0, "T":0, "R":0, "Y":0, "S":0, "W":0, "K":0, "M":0, "B":0, "D":0, "H":0, "V":0}
   for nucleotide in seq_list:
      counter[nucleotide] +=1
   return counter

def process_gc_from_pipe(in_queue, window_size):
   """
   Process fasta from a multiprocessing.queue and
   compute the percentage of GC within a given window
   """
   base_counter  = 1
   window_size   = float(window_size)
   while True:
      seq_list   = in_queue.get()
      #print seq_list
      if seq_list != 'EOF':
         if type(seq_list) == str:
            if seq_list.startswith('>'):
               chromosome = seq_list.strip('>').split(" ")[0]
               print "variableStep chrom=" + chromosome + " span=" + str(int(window_size))
               base_counter = 1
         else:
            seq_list = filter(None,seq_list)
            stats    = seq_map(seq_list)
            seq_len  = len(seq_list)
            if seq_len == window_size:
               if (stats['N'] + stats['R'] + stats['Y'] + stats['K'] + stats['M'] + stats['B'] + stats['D'] + stats['H'] + stats['V'])/window_size >= 0.75:
                  base_counter = base_counter + window_size
                  pass
               else:
                  gc_percent = 100 * ((stats['G'] + stats['C'] + stats['S'])/window_size)
                  print str(int(base_counter)) + "\t" + str(gc_percent)
                  base_counter = base_counter + window_size
            else:
               if (stats['N'] + stats['R'] + stats['Y'] + stats['K'] + stats['M'] + stats['B'] + stats['D'] + stats['H'] + stats['V'])/window_size >= 0.75:
                  base_counter = base_counter + seq_len
                  pass
               else:
                  gc_percent   = 100 * ((stats['G'] + stats['C'] + stats['S'])/float(seq_len))
                  print str(int(base_counter)) + "\t" + str(gc_percent)
                  base_counter = base_counter + seq_len
      else:
         break

class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _format_action(self, action):
        parts = super(argparse.RawDescriptionHelpFormatter, self)._format_action(action)
        if action.nargs == argparse.PARSER:
            parts = "\n".join(parts.split("\n")[1:])
        return parts

class DefaultHelpParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def DOpup2seqz(p1, p2, gc, n2, n, qlimit, qformat, hom, het, fileout, out_header):
   stream_mpileup = IterableQueue()
   if not n2:
      line_worker_partial = partial(line_worker, depth_sum=n, qlimit=qlimit, qformat=qformat, hom_t=hom, het_t=het)
      with xopen(p1, 'rb') as normal, xopen(p2, 'rb') as tumor, xopen(gc, 'rb') as gc_file:
         pup = multiPileups(normal,tumor)
         pup = GCmultiPileups(pup, gc_file)
         fileout.write("\t".join(out_header) + '\n')
         for line in pup:
            res = line_worker_partial(line)
            if res:
               fileout.write('\t'.join(map(str,res))+'\n')
   else:
      line_worker_partial = partial(line_worker, depth_sum=n, qlimit=qlimit, qformat=qformat, hom_t=hom, het_t=het, alt_pileup=True) 
      with xopen(p1, 'rb') as normal, xopen(p2, 'rb') as tumor, xopen(gc, 'rb') as gc_file, xopen(n2, 'rb') as alt_normal:
         pup = multiPileups(normal,tumor)
         pup = GCmultiPileups(pup, gc_file)
         pup = GCmultiPileupsAltDepth(pup, alt_normal)
         fileout.write("\t".join(out_header) + '\n')
         for line in pup:
            res = line_worker_partial(line)
            if res:
               fileout.write('\t'.join(map(str,res))+'\n')


def pileup2acgt(parser, subparser):
   subparser.add_argument('pileup', metavar='pileup',
                   help='Name of the input pileup (SAMtools) file. If the filename ends in .gz it will be opened in gzip mode. If the file name is - it will be read from STDIN.')
   subparser.add_argument('-n', dest='n', type=int,
                   help='The minimum read depth on a base required to make a mutation call.')
   parser_pup2muoutput      = subparser.add_argument_group(title='Output', description='Arguments that involve the output destination.')
   parser_pup2muoutput.add_argument('-o', '--output', dest='output', type = str, default = '-',
                       help='Name of the output file. To use gzip compression name the file ending in .gz. (default STDOUT).')
   parser_pup2muoutput.add_argument('--quiet', dest='quiet', action="store_true",
                       help='Do not output additional debugging information.')
   parser_pup2muqualitysets = subparser.add_argument_group(title='Quality and Format', description='Argument that change the quality threshold or the quality format.')
   parser_pup2muqualitysets.add_argument('-q', '--qlimit', dest='qlimit', default=20,type=int,
                   help='Minimum nucleotide quality score for inclusion in the counts.')
   parser_pup2muqualitysets.add_argument('-f', '--qformat', dest='qformat', default="sanger",
                   help='Quality format, options are "sanger" or "illumina". This will add an offset of 33 or 64 respectively to the qlimit value.')
   parser_pup2muinfo        = subparser.add_argument_group(title='Info',description='Other Information.')
   parser_pup2muinfo.add_argument('-v', '--version', action="version", version='%(prog)s version: ' + VERSION + ". " + DATE,
                   help='Display the version information and exit.')
   return parser.parse_args()

def pileup2seqz(parser, subparser):
   parser_ABinput    = subparser.add_argument_group(title='Input Files',description='Required input files.')
   parser_ABgenotype    = subparser.add_argument_group(title='Genotyper',description='Options regarding the genotyping.')
   parser_ABqualitysets = subparser.add_argument_group(title='Quality and Format', description='Options that change the quality threshold and format.')
   parser_ABinput.add_argument('-n', '--normal', dest = 'normal', required = True,
                   help='Name of the pileup file from the reference/normal sample')
   parser_ABinput.add_argument('-t', '--tumor', dest = 'tumor', required = True,
                   help='Name of the pileup file from the tumor sample')
   parser_ABinput.add_argument('-gc', dest = 'gc', metavar = 'gc', required = True,
                   help='The GC-content file coming from UCSC genome browser, or generated in the same UCSC format')
   parser_ABinput.add_argument('-n2', '--normal2', dest = 'normal2', type = str, default = None,
                   help='EXPERIMENTAL: Optional pileup used only to compute the depth.normal and depth-ratio, instead of \"normal\"')
   parser_ABqualitysets.add_argument('-q', '--qlimit', dest = 'qlimit', default = 20, type = int,
                   help='Minimum nucleotide quality score for inclusion in the counts. Default 20.')
   parser_ABqualitysets.add_argument('-f', '--qformat', dest = 'qformat', default = "sanger",
                   help='Quality format, options are "sanger" or "illumina". This will add an offset of 33 or 64 respectively to the qlimit value. Default "sanger".')
   parser_ABqualitysets.add_argument('-N', dest = 'n', type = int, default = 20,
                   help='Threshold to filter positions by the sum of read depth of the two samples. Default 20.')
   parser_ABgenotype.add_argument('--hom', dest = 'hom', type = float, default = 0.9,
                   help='Threshold to select homozygous positions. Default 0.9.')
   parser_ABgenotype.add_argument('--het', dest = 'het', type = float, default = 0.25,
                   help='Threshold to select heterozygous positions. Default 0.25.')
   return parser.parse_args()

def bam2seqz(parser, subparser):
   parser_ABinput       = subparser.add_argument_group(title='Input Files',description='Required input files.')
   parser_ABgenotype    = subparser.add_argument_group(title='Genotyper',description='Options regarding the genotyping.')
   parser_ABsamtools = subparser.add_argument_group(title='Samtools', description='Options regarding samtools.')
   parser_ABqualitysets = subparser.add_argument_group(title='Quality and Format', description='Options that change the quality threshold and format.')
   parser_ABinput.add_argument('-n', '--normal', dest = 'normal', required = True,
                   help='Name of the BAM file from the reference/normal sample')
   parser_ABinput.add_argument('-t', '--tumor', dest = 'tumor', required = True,
                   help='Name of the BAM file from the tumor sample')
   parser_ABinput.add_argument('-gc', dest = 'gc', metavar = 'gc', required = True,
                   help='The GC-content file coming from UCSC genome browser, or generated in the same UCSC format')
   parser_ABinput.add_argument('-F', "--fasta", dest = 'fasta', metavar = 'fasta', required = True,
                   help='The reference FASTA file used to generate the pileup')
   parser_ABinput.add_argument('-n2', '--normal2', dest = 'normal2', type = str, default = None,
                   help='EXPERIMENTAL: Optional BAM used only to compute the depth.normal and depth-ratio, instead of \"normal\"')
   parser_ABqualitysets.add_argument('-q', '--qlimit', dest = 'qlimit', default = 20, type = int,
                   help='Minimum nucleotide quality score for inclusion in the counts. Default 20.')
   parser_ABqualitysets.add_argument('-f', '--qformat', dest = 'qformat', default = "sanger",
                   help='Quality format, options are "sanger" or "illumina". This will add an offset of 33 or 64 respectively to the qlimit value. Default "sanger".')
   parser_ABqualitysets.add_argument('-N', dest = 'n', type = int, default = 20,
                   help='Threshold to filter positions by the sum of read depth of the two samples. Default 20.')
   parser_ABgenotype.add_argument('--hom', dest = 'hom', type = float, default = 0.9,
                   help='Threshold to select homozygous positions. Default 0.9.')
   parser_ABgenotype.add_argument('--het', dest = 'het', type = float, default = 0.25,
                   help='Threshold to select heterozygous positions. Default 0.25.')
   parser_ABsamtools.add_argument("-S", '--samtools', dest = 'samtools', type = str, default = "samtools",
                   help='Path of samtools to use for the pileup generation.')
   parser_ABsamtools.add_argument("-C", '--chromosome', dest = 'chr', type = str, default = None,
                   help='Argument to restrict the input/output to a chromosome or a chromosome region. Coordinate format is Name:pos.start-pos.end, eg: chr17:7565097-7590856, for a particular region; eg: chr17, for the entire chromosome. Chromosome names can checked in the BAM files and are depending on the FASTA reference used for alignment. Default behaviour is to not selecting any cromosome.')
   return parser.parse_args()

def GC_windows(parser, subparser):
   subparser.add_argument('fasta', metavar = 'fasta',
                   help='the fasta file. It can be a file name or \"-\" to take the input from STDIN')
   subparser.add_argument('-w', dest = 'window', type = int, default = 50,
                   help='The window size to calculate the GC-content percentage')
   return parser.parse_args()

def merge_pileups(parser, subparser):
   subparser.add_argument('-1', '--pileup1', dest = 'p1', required = True,
                   help='Name of the first pileup file')
   subparser.add_argument('-2', '--pileup2', dest = 'p2', required = True,
                   help='Name of the second pileup, will show as the last columns set')
   return parser.parse_args()

def reduce_seqz(parser, subparser):
   subparser.add_argument('-s', '--seqz', dest = 'seqz', required = True,
                   help='An seqz file from the pileup2seqz function.')
   subparser.add_argument('-w', '--window', dest = 'w', type = int, default = 50,
                   help='Window size used for binning the original seqz file. Default is 50.')
   return parser.parse_args()

def main():
   '''
   Execute the function with args
   '''
   parser = DefaultHelpParser(prog = __file__, formatter_class=lambda prog: SubcommandHelpFormatter(prog, max_help_position=20, width=75),
                              description='\nsequenza-utils.py \nThis script is part of Sequenza http://www.cbs.dtu.dk/biotools/sequenza/ \n\nSequenza Utils is a set of tools to perform various tasks, primarily aimed to convert bam/pileup files to a format usable by the sequenza R package.',
                              usage= 'sequenza-utils.py command [options]', epilog = 'This is version {0} - Francesco Favero - {1}'.format(VERSION, DATE))
   subparsers = parser.add_subparsers(dest='module')
   subparsers.metavar = None
   parser_bam2seqz  = subparsers.add_parser('bam2seqz', help = 'Process a paired set of BAM files (tumor and matching normal), and GC-content genome-wide information, to extract the common positions with A and B alleles frequencies',formatter_class=lambda prog: SubcommandHelpFormatter(prog,max_help_position=39, width=90))
   parser_pileup2seqz  = subparsers.add_parser('pileup2seqz', help = 'Process a paired set of pileup files (tumor and matching normal), and GC-content genome-wide information, to extract the common positions with A and B alleles frequencies',formatter_class=lambda prog: SubcommandHelpFormatter(prog,max_help_position=39, width=90))
   parser_reduce_seqz = subparsers.add_parser('seqz-binning', help = 'Perform binning of the seqz file to reduce file size and memory requirement for the analysis.')
   parser_pup2mu = subparsers.add_parser('pileup2acgt', help = 'Convert pileup format to ACGT format',formatter_class=lambda prog: SubcommandHelpFormatter(prog,max_help_position=30, width=90))
   parser_gc_window  = subparsers.add_parser('GC-windows', help = 'Compute the average GC content in sliding windows from a fasta file, output in the same format as gc5Base from UCSC')
   parser_merge_pileups = subparsers.add_parser('merge-pileups', help = 'Merge two pileups by finding the common positions, and return an mpileup file adding the second pileup as the last 3 columns.')
   try:
      used_module =  sys.argv[1]
      if used_module == "pileup2acgt":
         args = pileup2acgt(parser, parser_pup2mu)

         with xopen(args.output, "wb") as fileout, xopen(args.pileup, "rb") as f:
            fileout.write('chr' + "\t" + 'n_base' + "\t" + 'ref_base' + "\t" +  'read.depth' + "\t" + 'A' + "\t" + 'C' + "\t" + 'G' + "\t" + 'T' + "\t" + "strand" + '\n')
            parse_pileup_partial = partial(parse_pileup_str, min_depth=args.n, qlimit=args.qlimit, qformat=args.qformat)
            for line in f:
               try:
                  fileout.write(parse_pileup_partial(line) + '\n')
               except AttributeError:
                  pass                        

      elif used_module == "bam2seqz":
         args = bam2seqz(parser, parser_bam2seqz)
         with xopen('-', "wb") as fileout:
            out_header = ["chromosome", "position", "base.ref", "depth.normal", "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal", "GC.percent", "good.reads", "AB.normal", "AB.tumor", "tumor.strand"]
            if args.chr:
               cmd_tum1 = [args.samtools, "view", "-u", args.tumor, args.chr]
               cmd_nor1 = [args.samtools, "view", "-u", args.normal, args.chr]
            else:
               cmd_tum1 = [args.samtools, "view", "-u", args.tumor]
               cmd_nor1 = [args.samtools, "view", "-u", args.normal]

            cmd_tum2 = [args.samtools, "mpileup", "-f", args.fasta, "-q " + str(args.qlimit), "-"]
            cmd_nor2 = [args.samtools, "mpileup", "-f", args.fasta, "-q " + str(args.qlimit), "-"]
            if args.normal2:
               if args.chr:
                  cmd_nor3 = [args.samtools, "view", "-u", args.normal2, args.chr]
               else:
                  cmd_nor3 = [args.samtools, "view", "-u", args.normal2]
               cmd_nor4 = [args.samtools, "mpileup", "-f", args.fasta, "-q " + str(args.qlimit), "-"]
               tum1 = Popen(cmd_tum1, stdout = PIPE)
               nor1 = Popen(cmd_nor1, stdout = PIPE)
               nor2 = Popen(cmd_nor3, stdout = PIPE)
               with named_pipe() as tfifo, named_pipe() as nfifo, named_pipe() as n2fifo:
                  res = multiprocessing.Process(target = DOpup2seqz, args = (nfifo, tfifo, args.gc, n2fifo, args.n,  args.qlimit, args.qformat, args.hom, args.het, fileout, out_header))
                  res.start()
                  with open(nfifo, 'wb', 0) as normal, open(tfifo, 'wb', 0) as tumor, open(n2fifo, 'wb', 0) as normal2:
                     fifos = [Popen(cmd_tum2, stdin=tum1.stdout, stdout=tumor, stderr=PIPE), Popen(cmd_nor2, stdin=nor1.stdout, stdout=normal, stderr=PIPE), Popen(cmd_nor4, stdin=nor2.stdout, stdout=normal2, stderr=PIPE)]
                     tum1.stdout.close()
                     nor1.stdout.close()
                     nor2.stdout.close()
                  res.join()
                  for fifo in fifos:
                     fifo.terminate()
            else:
               tum1 = Popen(cmd_tum1, stdout = PIPE)
               nor1 = Popen(cmd_nor1, stdout = PIPE)
               with named_pipe() as tfifo, named_pipe() as nfifo:
                  res = multiprocessing.Process(target = DOpup2seqz, args = (nfifo, tfifo, args.gc, args.normal2, args.n,  args.qlimit, args.qformat, args.hom, args.het, fileout, out_header))
                  res.start()
                  with open(nfifo, 'wb', 0) as normal, open(tfifo, 'wb', 0) as tumor:
                     fifos = [Popen(cmd_tum2, stdin=tum1.stdout, stdout=tumor, stderr=PIPE), Popen(cmd_nor2, stdin=nor1.stdout, stdout=normal, stderr=PIPE)]
                     tum1.stdout.close()
                     nor1.stdout.close()
                  res.join()
                  for fifo in fifos:
                     fifo.terminate()

      elif used_module == "pileup2seqz":
         args = pileup2seqz(parser, parser_pileup2seqz)
         with xopen('-', "wb") as fileout:
            out_header = ["chromosome", "position", "base.ref", "depth.normal", "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal", "GC.percent", "good.reads", "AB.normal", "AB.tumor", "tumor.strand"]
            res = multiprocessing.Process(target = DOpup2seqz, args = (args.normal, args.tumor, args.gc, args.normal2, args.n,  args.qlimit, args.qformat, args.hom, args.het, fileout, out_header))
            res.start()
         res.join()


      elif used_module == "GC-windows":
         args = GC_windows(parser, parser_gc_window)
         gc_queue   = SimpleQueue()
         gc_process = multiprocessing.Process(target = process_gc_from_pipe, args = (gc_queue, args.window))
         gc_process.deamon = True
         gc_process.start()
         with xopen(args.fasta, 'rb') as fa_file:
            stream_fasta(fa_file, args.window, gc_queue)
         #gc_process.terminate()
         #gc_process.join()
      elif used_module == "merge-pileups":
         args = merge_pileups(parser, parser_merge_pileups)
         with xopen(args.p1, 'rb') as pileup1, xopen(args.p2, 'rb') as pileup2:
            pup = multiPileups(pileup1,pileup2)
            for line in pup:
               print('\t'.join(map(str, line)) + '\n')
      elif used_module == "seqz-binning":
         args = reduce_seqz(parser, parser_reduce_seqz)
         with xopen(args.seqz, 'rb') as seqz:
            abfred = abfreReduce(seqz, args.w)
            print abfred._header
            for a in abfred:
               if a:
                  print '\t'.join(a['top'])
                  for mid in a['middle']:
                     print '\t'.join(mid)
                  if ['end'] != ['top']:
                     print '\t'.join(a['end'])
      else:
         return parser.parse_args()

   except IndexError:
      args = parser.parse_args()

if __name__ == "__main__":
   main()
