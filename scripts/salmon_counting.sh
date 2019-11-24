#!/bin/bash

##############################################################################
# Options:
#   -i Indexed trascriptome folder (input)
#   -f FIRST INPUT FASTQ (trimmed FASTQ file)
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED) (trimmed FASTQ file)
# 	-t NUMBER OF THREADS
#		-o OUTPUT file
##############################################################################
while getopts ":i:f:s:t:o:" opt; do
  case $opt in
  i) INDEXED_FASTA=$OPTARG ;;
  f) INPUT_1=$OPTARG ;;
  s) INPUT_2=$OPTARG ;;
  t) THREADS=$OPTARG ;;
  o) OUTPUT=$OPTARG ;;
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
if [ -z $THREADS ]; then
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

# Check indexed fasta file
if [ ! -d "$INDEXED_FASTA" ]; then
  echo "Indexed trascriptome does not exist!"
  exit 7
fi

#### Counting ####
OUTPUT_DIR=$(dirname "$OUTPUT")
TEMP_DIR="$OUTPUT_DIR/TMP"

if [ $PAIRED = "true" ]; then
  salmon quant -i "$INDEXED_FASTA" -l A -1 "$INPUT_1" -2 "$INPUT_2" --validateMappings -p $THREADS -o "$TEMP_DIR"
else
  salmon quant -i "$INDEXED_FASTA" -l A -r "$INPUT_1" --validateMappings -p $THREADS -o "$TEMP_DIR"
fi

OUTPUT_FILE="$TEMP_DIR/quant.sf"

# Check output file
if [ ! -f "$OUTPUT_FILE" ]; then
  echo "Unable to find output file!" >&2
  exit 8
fi

# Move output file from tmp directory to output directory
mv "$OUTPUT_FILE" "$OUTPUT"
chmod 777 "$OUTPUT"

# Removing items of tmp directory
if [ -d "$TEMP_DIR" ]; then
  rm -rf "$TEMP_DIR"
fi
