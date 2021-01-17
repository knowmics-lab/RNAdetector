#!/bin/bash

REFERENCES=$(find "/rnadetector/ws/storage/app/references" -name '*.fa' -or -name '*.fasta')

for FASTA in $REFERENCES; do
  if [ ! -f "$FASTA.fai" ]; then
    echo "No fasta index found for $FASTA. Indexing..."
    samtools faidx "$FASTA"
  fi
done

ANNOTATIONS=$(find "/rnadetector/ws/storage/app/annotations" -name '*.gtf')

for GTF in $ANNOTATIONS; do
  OUTPUT_DIRECTORY=$(dirname "$GTF")
  NAME=$(basename "$GTF" ".gtf")
  if [ ! -f "$OUTPUT_DIRECTORY/$NAME.gff3.gz" ]; then
    echo "No index found for $GTF. Indexing..."
    bash /rnadetector/scripts/prepare_gtf.sh "$GTF"
  fi
done
