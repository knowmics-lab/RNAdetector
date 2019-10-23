#!/bin/bash

# Path directories
FASTQ_PATH=  # Path della cartella con i file FASTQ di input
QUALITY_FILTERED_FASTQ_PATH=  # Path della cartella con i file FASTQ in cui sono state filtrati le reads di bassa qualità e tolti gli adattatori
FASTQC_REPORTS_PATH=  # Path della cartella contenente i reports di FASTQC
CONVERTED_FASTQ_PATH=  # Path dell cartella contenente i file FASTQ convertiti da file BAM
BAM_PATH=  # Path della cartella con i file BAM di input
ALIGNMENT_PATH= # Path della cartella contenete i file BAM in seguito all'allineamento
READ_COUNTS_PATH=  # Path della cartella che contiene i file di testo con le conte grezze

# Check directories
if [ ! -d $FASTQ_PATH ]; then
	mkdir -p $FASTQ_PATH
fi

if [ ! -d $QUALITY_FILTERED_FASTQ_PATH ]; then
	mkdir -p $QUALITY_FILTERED_FASTQ_PATH
fi

if [ ! -d $FASTQC_REPORTS_PATH ]; then
	mkdir -p $FASTQC_REPORTS_PATH
fi

if [ ! -d $CONVERTED_FASTQ_PATH ]; then
	mkdir -p $CONVERTED_FASTQ_PATH
fi

if [ ! -d $BAM_PATH ]; then
	mkdir -p $BAM_PATH
fi

if [ ! -d $ALIGNMENT_PATH ]; then
	mkdir -p $ALIGNMENT_PATH
fi

if [ ! -d $READ_COUNTS_PATH ]; then
	mkdir -p $READ_COUNTS_PATH
fi
