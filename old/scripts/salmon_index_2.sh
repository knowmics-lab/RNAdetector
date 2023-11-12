#!/bin/bash

##############################################################################
# Options:
# 	-r FASTA transcripts files
#   -i Indexed trascriptome (output folder)
##############################################################################
OTHER_ARGS=""
while getopts ":r:i:A:" opt; do
  case $opt in
  r) FASTA=$OPTARG ;;
  i) INDEXED_FASTA=$OPTARG ;;
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
# Check FASTA file
if [ -z "$FASTA" ] || [ ! -f "$FASTA" ]; then
  echo "FASTA file with transcripts does not exist!"
  exit 3
fi

[ ! -f "$FASTA.fai" ] && samtools faidx "$FASTA" && chmod 777 "$FASTA.fai"

if [ -n "$OTHER_ARGS" ]; then
  echo "Processing with custom arguments: \"${OTHER_ARGS}\""
fi

#### Indexed transcriptome ####
# shellcheck disable=SC2086
if ! salmon index -t "$FASTA" -i "$INDEXED_FASTA" -k 31 $OTHER_ARGS 2>/dev/null; then
  echo "An error occurred during salmon index execution!"
  exit 5
fi

if [ ! -d "$INDEXED_FASTA" ]; then
  echo "Indexed trascriptome does not exist!"
  exit 4
fi

chmod -R 777 "$INDEXED_FASTA"
