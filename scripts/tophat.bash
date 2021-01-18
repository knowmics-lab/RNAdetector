#!/bin/bash

##############################################################################
# Options:
# 	-a FILE GTF
# 	-g REFERENCE GENOME FILE (bowtie indexed genome path with prefix)
# 	-t NUMBER OF THREADS
# 	-f FIRST INPUT FASTQ (trimmed FASTQ file)
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED) (trimmed FASTQ file)
# 	-o OUTPUT BAM FILE
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
TEMP_DIR=$(mktemp -d)
if [ $PAIRED = "true" ]; then
  if ! tophat2 -G "$GTF_FILE" -o "$TEMP_DIR" -p "$THREADS" "$REF_GENOME" "$INPUT_1" "$INPUT_2"; then
    echo "An error occurred during tophat2 execution!"
    exit 8
  fi
else
  if ! tophat2 -G "$GTF_FILE" -o "$TEMP_DIR" -p "$THREADS" "$REF_GENOME" "$INPUT_1"; then
    echo "An error occurred during tophat2 execution!"
    exit 8
  fi
fi

BAM="$TEMP_DIR/accepted_hits.bam"

# Check BAM file
if [ ! -f "$BAM" ]; then
  echo "Unable to find output bam file!"
  exit 9
fi

bash /rnadetector/scripts/prepare_bam.sh -f "$OUTPUT" -s -u "$BAM" -t "$THREADS" || exit_abnormal "Unable to prepare BAM file" "$?"

# Removing items of tmp directory
if [ -d "$TEMP_DIR" ]; then
  rm -rf "$TEMP_DIR"
fi
