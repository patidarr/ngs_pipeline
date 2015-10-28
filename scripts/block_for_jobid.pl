#!/usr/bin/perl

# This script blocks until the jobid given is either unknown or not running/pending

use lib "/usr/local/slurm/lib/site_perl/5.12.1/x86_64-linux-thread-multi";
use lib "/usr/local/slurm/lib/perl5/site_perl/5.18.2/x86_64-linux-thread-multi-ld";
use Slurm;

my $known;
my $done;
sleep 5;
while (!$done) {
  my $slurm = Slurm::new();
  $jobs = $slurm->load_jobs();
  JOB: foreach my $ref (@{$jobs->{job_array}}) {
    next JOB unless ($ref->{job_id} == $ARGV[0]);
    if (defined $ref->{job_state}) { $known = 1; }
    if ($ref->{job_state} > 1)  { print "jobid $ARGV[0] ended\n"; $done=1; }
  }
  if (!$known) { print "jobid $ARGV[0] not known\n";  exit; }
  sleep 5;
}
