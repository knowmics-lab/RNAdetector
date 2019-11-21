#!/usr/bin/env bash

CURR=$(pwd)
git clone https://github.com/alessandrolaferlita/RNAdetector.git
tar -zcvf repo.tar.gz ./RNAdetector
rm -rf ./RNAdetector
mkdir tmp
cd tmp
wget --timestamping "https://alpha.dmi.unict.it/~alaimos/RNAdetector/hg19.fa.gz" -O "hg19.fa.gz"
wget --timestamping "https://alpha.dmi.unict.it/~alaimos/RNAdetector/hg19_circRNA_circbase.gtf.gz" -O "hg19_circRNA_circbase.gtf.gz"
wget --timestamping "https://alpha.dmi.unict.it/~alaimos/RNAdetector/hg19_small_ncRNA.gtf.gz" -O "hg19_small_ncRNA.gtf.gz"
gunzip -v *.gz
mv hg19.fa "$CURR/Human_hg19_genome.fasta"
mv hg19_circRNA_circbase.gtf "$CURR/Human_hg19_circRNAs.gtf"
mv hg19_small_ncRNA.gtf "$CURR/Human_hg19_small_ncRNAs.gtf"
cd $CURR
rm -r tmp
docker build -t alaimos/ubuntu-base .
rm repo.tar.gz
rm *.gtf
rm *.fasta