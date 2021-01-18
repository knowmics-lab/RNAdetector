#!/bin/bash

##############################################################################
# Options:
# 	-f FASTA input genome file
#   -p STR	The directory containing the index
# 	-t NUMBER OF THREADS
##############################################################################
while getopts ":f:p:t:" opt; do
  case $opt in
  f) FASTA_FILE=$OPTARG ;;
  p) PREFIX_OUTPUT=$OPTARG ;;
  t) THREADS=$OPTARG ;;
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
# Check number of threads and set 1 as default value
if [ -z "$THREADS" ]; then
  THREADS=1
fi

# Check input files
if [ -z "$FASTA_FILE" ] || [ ! -f "$FASTA_FILE" ]; then
  echo "Input file does not exist!"
  exit 3
fi

# Check output
if [ -z "$PREFIX_OUTPUT" ]; then
  echo "Output directory must be specified!"
  exit 4
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$PREFIX_OUTPUT")" ]; then
  echo "Output directory is not writable!"
  exit 5
fi

[ ! -f "$FASTA_FILE.fai" ] && samtools faidx "$FASTA_FILE" && chmod 777 "$FASTA_FILE.fai"

REFERENCE_DIR="${PREFIX_OUTPUT}_star"

[ ! -d "$REFERENCE_DIR" ] && mkdir "$REFERENCE_DIR"

#### Genome indexing ####
STAR --runThreadN "$THREADS" \
  --runMode genomeGenerate \
  --genomeDir "$REFERENCE_DIR" \
  --genomeFastaFiles "$FASTA_FILE" || (
  echo "An error occurred during indexing!" && exit 6
)

chmod -R 777 "$REFERENCE_DIR"
