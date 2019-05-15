#!/bin/bash

##############################################################################
# Options:
#		-a ADAPTER (OPZIONALE) --illumina --nextera --small_rna
#		-q QUALITY (default 20)
#		-l LENGHT (default 14)
# 	-f FIRST INPUT FASTQ
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED)
# 	-o OUTPUT TRIMMED FASTQ
##############################################################################
while getopts ":a:q:l:f:s:o:" opt; do
	case $opt in
		a ) ADAPTER=$OPTARG ;;
		q ) QUALITY=$OPTARG ;;
		l ) LENGHT=$OPTARG ;;
		f ) INPUT_1=$OPTARG ;;
		s ) INPUT_2=$OPTARG ;;
		o ) OUTPUT=$OPTARG ;;
		\?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
		: ) echo "Option -$OPTARG requires an argument." >&2; exit 2;;
	esac
done

#### Check parameters ####
# Quality control
if [ -z $QUALITY ]; then
	$QUALITY = 20
fi
# Control read lenght
if [ -z $LENGHT ]; then
	$LENGHT = 14
fi
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

#### Trimming and adaptors removing ####
if [ $PAIRED ]; then
	trim_galore -q $QUALITY --paired -a $ADAPTER -a2 $ADAPTER -o $OUTPUT --dont_gzip --phred33 --length $LENGHT --no_report_file  $INPUT_1 $INPUT_2
else
	trim_galore -q $QUALITY -a $ADAPTER -o $OUTPUT --dont_gzip --phred33 --length $LENGHT --no_report_file $INPUT_1
fi
