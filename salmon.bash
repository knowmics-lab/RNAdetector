#!/bin/bash

##############################################################################
# Options:
# 	-r FASTA transcripts files
#   -i Indexed trascriptome
#   -d decoys file
#   -f FIRST INPUT FASTQ (trimmed FASTQ file)
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED) (trimmed FASTQ file)
# 	-t NUMBER OF THREADS
#		-o OUTPUT
##############################################################################
while getopts ":r:i:d:f:s:t:o:" opt; do
	case $opt in
		r ) FASTA=$OPTARG ;;
		i ) INDEXED_FASTA=$OPTARG ;;
    d ) DECOYS_FILE=$OPTARG ;;
		f ) INPUT_1=$OPTARG ;;
		s ) INPUT_2=$OPTARG ;;
		t ) THREADS=$OPTARG ;;
		o ) OUTPUT=$OPTARG ;;
		\?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
		: ) echo "Option -$OPTARG requires an argument." >&2; exit 2;;
	esac
done

#### Check parameters ####

# Check FASTA file
if [ -z $FASTA ] || [ ! -f $FASTA ]; then
	echo "FASTA file does not exist!"
	exit 3
fi

# Check input files
if [ -z $INPUT_1 ] || [ ! -f $INPUT_1 ]; then
	echo "Input file does not exist!" >&2
	exit 4
fi

# Control sequencing strategy "single end" o "paired end"
if [ -z $INPUT_2 ]; then
	PAIRED=false
elif [ ! -f $INPUT_2 ]; then
	echo "Second input file does not exist!" >&2
	exit 5
else
    PAIRED=true
fi

# Check number of threads and set 1 as default value
if [ -z $THREADS ]; then
	THREADS=1
fi

# Check output
if [ -z $OUTPUT ]; then
	echo "Output file must be specified!" >&2
	exit 6
fi

# Check if output directory is writable
if [ ! -w $(dirname $OUTPUT) ]; then
	echo "Output directory is not writable!" >&2
	exit 7
fi

# Check if decoy file exists
# Check FASTA file
if [ -z $DECOYS_FILE ] || [ ! -f $DECOYS_FILE ]; then
	echo "Decoys file does not exist!"
	exit 8
fi

#### Indexed transcriptome ####
if [ -z $INDEXED_FASTA ]; then
	sudo docker run combinelab/salmon salmon index -t $FASTA -i $INDEXED_FASTA -decoys $DECOYS_FILE -k 31
fi

# Check indexed fasta file
if [ -z $INDEXED_FASTA ] || [ ! -f $INDEXED_FASTA ]; then
	echo "Indexed FASTA file does not exist!"
	exit 9
fi

#### Counting ####
if [ $PAIRED ]; then
	sudo docker run combinelab/salmon salmon quant -i $INDEXED_FASTA -l A -1 $INPUT_1 -2 $INPUT_2 --validateMappings -p $THREADS -o $OUTPUT
else
  sudo docker run combinelab/salmon salmon quant -i $INDEXED_FASTA -l A -r $INPUT_1 --validateMappings -p $THREADS -o $OUTPUT
fi
