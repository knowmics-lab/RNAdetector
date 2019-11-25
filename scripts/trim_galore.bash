#!/bin/bash

##############################################################################
# Options:
#		-q QUALITY (default 20)
#		-l LENGHT (default 14)
# 	-f FIRST INPUT FASTQ
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED)
# 	-t NUMBER OF THREADS
# 	-o OUTPUT TRIMMED FASTQ
##############################################################################

HARD_TRIM=false

while getopts ":q:l:f:s:t:o:h" opt; do
  case $opt in
  q) QUALITY=$OPTARG ;;
  l) LENGHT=$OPTARG ;;
  f) INPUT_1=$OPTARG ;;
  s) INPUT_2=$OPTARG ;;
  t) THREADS=$OPTARG ;;
  o) OUTPUT=$OPTARG ;;
  h) HARD_TRIM=true ;;
  \?)
    echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
  :)
    echo "Option -$OPTARG requires an argument." >&2
    exit 2
    ;;
  esac
done

#### Check parameters ####
# Quality control
if [ -z "$QUALITY" ]; then
  QUALITY=20
fi
# Control read lenght
if [ -z "$LENGHT" ]; then
  LENGHT=14
fi
# Check input files
if [ -z "$INPUT_1" ] || [ ! -f "$INPUT_1" ]; then
  echo "Input file does not exist!" >&2
  exit 3
fi
# Control sequencing strategy "single end" o "paired end"
if [ -z "$INPUT_2" ]; then
  PAIRED=false
elif [ ! -f "$INPUT_2" ]; then
  echo "Second input file does not exist!" >&2
  exit 4
else
  PAIRED=true
fi

# Check number of threads and set 1 as default value
if [ -z "$THREADS" ]; then
  THREADS=1
fi

# Check output
if [ -z "$OUTPUT" ]; then
  echo "Output file must be specified!" >&2
  exit 5
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$OUTPUT")" ]; then
  echo "Output directory is not writable!" >&2
  exit 6
fi

#### Trimming and adaptors removing ####
if [ $PAIRED = "true" ]; then
  trim_galore -j $THREADS -q $QUALITY --paired -o "$OUTPUT" --dont_gzip --phred33 --length $LENGHT --no_report_file "$INPUT_1" "$INPUT_2"
  if [ $HARD_TRIM = "true" ]; then
    trim_galore -j $THREADS -q $QUALITY --paired -o "$OUTPUT" --dont_gzip --phred33 --length $LENGHT --hardtrim5 $LENGHT --no_report_file "$INPUT_1" "$INPUT_2"
    #TODO
  fi
else
  trim_galore -j $THREADS -q $QUALITY -o "$OUTPUT" --dont_gzip --phred33 --length $LENGHT --no_report_file "$INPUT_1"
  if [ $HARD_TRIM = "true" ]; then
    trim_galore -j $THREADS -q $QUALITY -o "$OUTPUT" --dont_gzip --phred33 --length $LENGHT --hardtrim5 $LENGHT --no_report_file "$INPUT_1"
    #TODO
  fi
fi
