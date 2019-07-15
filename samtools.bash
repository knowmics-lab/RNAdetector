#!/bin/bash

##############################################################################
# Options:
# 	-b INPUT BAM
# 	-o OUTPUT FASTQ
##############################################################################
while getopts ":b:o:" opt; do
	case $opt in
		b ) INPUT=$OPTARG ;;
		o ) OUTPUT=$OPTARG ;;
		\?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
		: ) echo "Option -$OPTARG requires an argument." >&2; exit 2;;
	esac
done

#### Check parameters ####
# Check input files
if [ -z $INPUT ] || [ ! -f $INPUT ]; then
	echo "Input file does not exist!" >&2
	exit 3
fi
# Check output
if [ -z $OUTPUT ]; then
	echo "Output file must be specified!" >&2
	exit 4
fi
# Check if output directory is writable
if [ ! -w $(dirname $OUTPUT) ]; then
	echo "Output directory is not writable!" >&2
	exit 5
fi

#### Conversion from BAM to FASTQ ####
samtools fastq $INPUT > $OUTPUT
