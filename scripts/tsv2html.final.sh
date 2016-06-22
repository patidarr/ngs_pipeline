#!/bin/bash
 
usage()
{
cat <<EOF
 
Usage: $(baselocation $0) [OPTIONS] input > output
 
Script to produce HTML tables from delimited input. Delimiter can be specified
as an optional argument. If omitted, script defaults to comma.
 
Options:
 
  -d       Specify delimiter to look for, instead of "\t".
 
  --head   Treat first line as header, enclosing in <thead> and <th> tags.
 
  --foot   Treat last line as footer, enclosing in <tfoot> and <th> tags. 

  --Name   To include in Mail
  --Diagnosis To include in the mail body  
Examples:
 
  1. $(baselocation $0) --Name NCI00001 input.csv
 
  Above will parse file 'input.csv' with comma as the field separator and
  output HTML tables to STDOUT.
 
  2. $(baselocation $0) -d '|' < input.psv > output.htm
 
  Above will parse file "input.psv", looking for the pipe character as the
  delimiter, then output results to "output.htm".
 
  3. $(baselocation $0) -d '\t' --head --foot < input.tsv > output.htm
 
  Above will parse file "input.tsv", looking for tab as the delimiter, then
  process first and last lines as header/footer (that contain data labels), then
  write output to "output.htm".
 
EOF
}
 
while true; do
  case "$1" in
    -d)
      shift
      d="$1"
      ;;
    --foot)
      foot="-v ftr=1"
      ;;
    --location)
      shift
      location="$1"
      ;;
    --host)
      shift
      host="$1"
      ;;
    --help)
      usage
      exit 0
      ;;
    --head)
      head="-v hdr=1"
      ;;
    -*)
      echo "ERROR: unknown option '$1'"
      echo "see '--help' for usage"
      exit 1
      ;;
    *)
      f=$1
      break
      ;;
  esac
  shift
done
 
if [ -z "$d" ]; then
  d="\t"
fi
 
if [ -z "$f" ]; then
  echo "ERROR: input file is required"
  echo "see '--help' for usage"
  exit 1
fi
 
if ! [ -f "$f" ]; then
  echo "ERROR: input file '$f' is not readable"
  exit 1
else
  data=$(sed '/^$/d' $f)
  last=$(wc -l <<< "$data")
fi
echo "<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<style type="text/css">
#main {
	padding: 5px;
	border-collapse: collapse; 
	border: 1px solid #000000;
	width: 100%;
}
#main td {
	border: 1px solid #000000;
	padding: 3px;
	font-size: .9em;
}
#main th {
	border: 1px solid #000000;
	#height: 100px;
	background-color: #CCFFCC;
	#transform:rotate(45deg);
}
</style>

</head>
<body>
<p>
Hello,<br><br>
      ngs-pipeline finished successfully on <b>$host</b><br><br> 
	

Subject(s) Processed:<br><br>
</p>
"
 
awk -F "$d" -v last=$last $head $foot '
  BEGIN {
    print "  <table id=\"main\">"
  }       
  {
    gsub(/</, "\\&lt;")
    gsub(/>/, "\\&gt;")
    if(NR == 1 && hdr) {  
      printf "    <thead>\n"
    gsub(/&/, "\\&gt;")    }
    if(NR == last && ftr) {  
      printf "    <tfoot>\n"
    }
    print "      <tr>"
    for(f = 1; f <= NF; f++)  {
      if((NR == 1 && hdr) || (NR == last && ftr)) {
        printf "        <th align=\"center\" >%s</th>\n", $f
      }
      else
		printf "        <td align=\"center\" >%s</td>\n", $f
    }     
    print "      </tr>"
    if(NR == 1 && hdr) {
      printf "    </thead>\n"
    }
    if(NR == last && ftr) {
      printf "    </tfoot>\n"
    }
  }       
  END {
    print "  </table>"
  }
' <<< "$data"
echo "
<br><br>
Result available in <b>$location</b>.
<br>
<p>For IGV sessions to work in igv:<br>
biowulf.nih.gov(/data/khanlab) is mapped as K <br>
biowulf.nih.gov(/data/Clinomics) is mapped as V <br>
login01 aks TGen Server(/projects/Clinomics) is mapped as Y<br><br>

If you ran pipeline on another location, igv session file can not be loaded in IGV.</p><br><br>

<p>Regards,<br>
KhanLab<br>
Oncogenomics Section<br>
CCR NCI NIH<br>
</p>
</body>
</html>"
