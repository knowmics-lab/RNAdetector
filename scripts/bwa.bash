#!/bin/bash

##############################################################################
# Options:
# 	-a FILE GTF
# 	-g BWA INDEXED REFERENCE GENOME FOLDER WITH PREFIX
# 	-t NUMBER OF THREADS
# 	-f FIRST INPUT FASTQ (trimmed FASTQ file)
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED) (trimmed FASTQ file)
# 	-o OUTPUT SAM FILE
##############################################################################
while getopts ":a:g:t:f:s:o:" opt; do
  case $opt in
  a) GTF_FILE=$OPTARG ;;
  g) REF_GENOME=$OPTARG ;;
  t) THREADS=$OPTARG ;;
  f) INPUT_1=$OPTARG ;;
  s) INPUT_2=$OPTARG ;;
  o) OUTPUT=$OPTARG ;;
  \?)
    echo "Invalid option: -$OPTARG"
    exit 1
    ;;
  :)
    echo "Option -$OPTARG requires an argument."
    exit 2
    ;;
  esac
done

#### Check parameters ####

# Check GTF annotation files
if [ -z "$GTF_FILE" ] || [ ! -f "$GTF_FILE" ]; then
  echo "Annotation file does not exist!"
  exit 3
fi

# Check input files
if [ -z "$INPUT_1" ] || [ ! -f "$INPUT_1" ]; then
  echo "Input file does not exist!"
  exit 4
fi

# Control sequencing strategy "single end" o "paired end"
if [ -z "$INPUT_2" ]; then
  PAIRED=false
elif [ ! -f "$INPUT_2" ]; then
  echo "Second input file does not exist!"
  exit 5
else
  PAIRED=true
fi

# Check number of threads and set 1 as default value
if [ -z "$THREADS" ]; then
  THREADS=1
fi

# Check output
if [ -z "$OUTPUT" ]; then
  echo "Output file must be specified!"
  exit 6
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$OUTPUT")" ]; then
  echo "Output directory is not writable!"
  exit 7
fi

#### Alignment ####
if [ $PAIRED = "true" ]; then
  if ! fastq_pair "$INPUT_1" "$INPUT_2"; then
    echo "Unable to determine unpaired reads"
    exit 9
  fi
  if [ ! -f "$INPUT_1.paired.fq" ]; then
    echo "Unable to find paired reads of first input file!"
    exit 10
  fi
  if [ ! -f "$INPUT_2.paired.fq" ]; then
    echo "Unable to find paired reads of second input file!"
    exit 11
  fi
  rm "$INPUT_1.single.fq"
  rm "$INPUT_2.single.fq"
  if ! bwa mem -t "$THREADS" "$REF_GENOME" "$INPUT_1.paired.fq" "$INPUT_2.paired.fq" >"$OUTPUT"; then
    echo "An error occurred during bwa mem execution!"
    exit 12
  fi
else
  if ! bwa mem -t "$THREADS" "$REF_GENOME" "$INPUT_1" >"$OUTPUT"; then
    echo "An error occurred during bwa mem execution!"
    exit 12
  fi
fi

#Check SAM file
if [ ! -f "$OUTPUT" ]; then
  echo "Unable to find bwa output file!"
  exit 8
fi

chmod 777 "$OUTPUT"
