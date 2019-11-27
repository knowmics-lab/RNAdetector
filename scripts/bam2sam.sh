#!/bin/bash

##############################################################################
# Options:
# 	-b INPUT BAM
# 	-o OUTPUT SAM FILE
##############################################################################
while getopts ":b:o:" opt; do
  case $opt in
  b) INPUT=$OPTARG ;;
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
# Check input files
if [ -z "$INPUT" ] || [ ! -f "$INPUT" ]; then
  echo "Input file does not exist!"
  exit 3
fi

# Check output
if [ -z "$OUTPUT" ]; then
  echo "Output file must be specified!"
  exit 4
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$OUTPUT")" ]; then
  echo "Output directory is not writable!"
  exit 5
fi

#### Conversione BAM to SAM ####
if ! samtools view -h "$INPUT" >"$OUTPUT"; then
  echo "An error occurred during samtools execution!"
  exit 6
fi
chmod 777 "$OUTPUT"
