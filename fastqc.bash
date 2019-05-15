#!/bin/bash

##############################################################################
# Options:
# 	-f FIRST INPUT FASTQ (trimmed FASTQ file)
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED) (trimmed FASTQ file)
# 	-o OUTPUT QUALITY REPORTS
##############################################################################
while getopts ":f:s:o:" opt; do
	case $opt in
		f ) INPUT_1=$OPTARG ;;
		s ) INPUT_2=$OPTARG ;;
		o ) OUTPUT=$OPTARG ;;
		\?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
		: ) echo "Option -$OPTARG requires an argument." >&2; exit 2;;
	esac
done

#### Check parameters ####
# Check input files
if [ -z $INPUT_1 ] || [ ! -f $INPUT_1 ]; then
	echo "Input file does not exist!" >&2
	exit 3
fi
# Control sequencing strategy "single end" o "paired end"
if [ -z $INPUT_2 ]; then
 PAIRED=false
elif [ ! -f $INPUT_2 ]; then
	echo "Second input file does not exist!" >&2
	exit 4
else
 PAIRED=true
fi
# Check output
if [ -z $OUTPUT ]; then
	echo "Output file must be specified!" >&2
	exit 5
fi
# Check if output directory is writable
if [ ! -w $(dirname $OUTPUT) ]; then
	echo "Output directory is not writable!" >&2
	exit 6
fi

#### Quality control with FASTQC ####
if [ $PAIRED ]; then
	fastqc $INPUT_1 $INPUT_2 --outdir=$OUTPUT
else
	fastqc $INPUT_1 --outdir=$OUTPUT
fi
