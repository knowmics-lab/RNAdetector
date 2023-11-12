#!/bin/bash

##############################################################################
# Options:
# 	-g REFERENCE INDEXED GENOME FILE (basename of indexed genome)
# 	-t NUMBER OF THREADS
# 	-f FIRST INPUT FASTQ (trimmed FASTQ file)
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED) (trimmed FASTQ file)
# 	-o OUTPUT BAM FILE
##############################################################################
OTHER_ARGS=""
while getopts ":g:t:f:s:o:A:" opt; do
  case $opt in
  g) REF_GENOME=$OPTARG ;;
  t) THREADS=$OPTARG ;;
  f) INPUT_1=$OPTARG ;;
  s) INPUT_2=$OPTARG ;;
  o) OUTPUT=$OPTARG ;;
  A) OTHER_ARGS=$OPTARG ;;
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
#Check input files
if [ -z "$INPUT_1" ] || [ ! -f "$INPUT_1" ]; then
  echo "Input file does not exist!"
  exit 3
fi

# Control sequencing strategy "single end" o "paired end"
if [ -z "$INPUT_2" ]; then
  PAIRED=false
elif [ ! -f "$INPUT_2" ]; then
  echo "Second input file does not exist!"
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
  echo "Output file must be specified!"
  exit 5
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$OUTPUT")" ]; then
  echo "Output directory is not writable!"
  exit 6
fi

if [ -n "$OTHER_ARGS" ]; then
  echo "Processing with custom arguments: \"${OTHER_ARGS}\""
fi

#### Alignment ####
TMP_OUTPUT="$(dirname "$OUTPUT")/temp.bam"
if [ $PAIRED = "true" ]; then
  # shellcheck disable=SC2086
  if ! hisat2 $OTHER_ARGS -p "$THREADS" -x "$REF_GENOME" -1 "$INPUT_1" -2 "$INPUT_2" | samtools view -Sbh >"$TMP_OUTPUT"; then
    echo "An error occurred during HISAT 2 execution!"
    exit 7
  fi
else
  # shellcheck disable=SC2086
  if ! hisat2 $OTHER_ARGS -p "$THREADS" -x "$REF_GENOME" -U "$INPUT_1" | samtools view -Sbh >"$TMP_OUTPUT"; then
    echo "An error occurred during HISAT 2 execution!"
    exit 7
  fi
fi

if [ ! -f "$TMP_OUTPUT" ]; then
  echo "Unable to find HISAT2 output file!"
  exit 8
fi

if ! bash /rnadetector/scripts/prepare_bam.sh -f "$OUTPUT" -s -u "$TMP_OUTPUT" -t "$THREADS"; then
  echo "Unable to find HISAT2 output file!"
  exit 9
fi

# Check SAM file
if [ ! -f "$OUTPUT" ]; then
  echo "Unable to find sorted output file!"
  exit 10
fi

rm "$TMP_OUTPUT"
