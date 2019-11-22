#!/usr/bin/env bash

service nginx start
service php7.2-fpm start
service supervisor start

if [ ! -f "/rnadetector/ws/storage/app/references/indexed" ]; then
    echo "Indexing genomes...this might take a while..."
    /bin/bash "/rnadetector/scripts/bwa_index.sh"     -f "/rnadetector/ws/storage/app/references/Human_hg19_genome/reference.fasta" -p "/rnadetector/ws/storage/app/references/Human_hg19_genome/reference"
    /bin/bash "/rnadetector/scripts/bowtie2_index.sh" -f "/rnadetector/ws/storage/app/references/Human_hg19_genome/reference.fasta" -p "/rnadetector/ws/storage/app/references/Human_hg19_genome/reference"
    /bin/bash "/rnadetector/scripts/salmon_index_2.sh" -r "/rnadetector/ws/storage/app/references/Human_hg19_transcriptome/reference.fasta" -i "/rnadetector/ws/storage/app/references/Human_hg19_transcriptome/reference"
    /bin/bash "/rnadetector/scripts/salmon_index_2.sh" -r "/rnadetector/ws/storage/app/references/Human_hg19_mRNA_transcriptome/reference.fasta" -i "/rnadetector/ws/storage/app/references/Human_hg19_mRNA_transcriptome/reference"
    /bin/bash "/rnadetector/scripts/salmon_index_2.sh" -r "/rnadetector/ws/storage/app/references/Human_hg19_lncRNA_transcriptome/reference.fasta" -i "/rnadetector/ws/storage/app/references/Human_hg19_lncRNA_transcriptome/reference"
    touch /rnadetector/ws/storage/app/references/indexed
fi

exec "$@"
