#!/bin/bash

REFERENCES=$(find "/rnadetector/ws/storage/app/references" -name '*.fa' -or -name '*.fasta')

for FASTA in $REFERENCES; do
  if [ ! -f "$FASTA.fai" ]; then
    echo "No fasta index found for $FASTA. Indexing..."
    samtools faidx "$FASTA"
  fi
  rm "$FASTA.fai" #@todo remove this line for production
done
