#!/bin/bash

##############################################################################
# Options:
#   -r FASTA transcripts files
#   -i input BAM file
# 	-t NUMBER OF THREADS
#		-o OUTPUT directory
##############################################################################
while getopts ":r:i:t:o:" opt; do
	case $opt in
	  r ) FASTA_TRANSCRIPTS=$OPTARG ;;
		i ) INPUT_BAM=$OPTARG ;;
    t ) THREADS=$OPTARG ;;
		o ) OUTPUT=$OPTARG ;;
		\?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
		: ) echo "Option -$OPTARG requires an argument." >&2; exit 2;;
	esac
done

#### Check parameters ####
# Check input files
if [ -z "$INPUT_BAM" ] || [ ! -f "$INPUT_BAM" ]; then
	echo "Input file does not exist!" >&2
	exit 3
fi

# Check check fasta transcriptome
if [ -z "$FASTA_TRANSCRIPTS" ] || [ ! -f "$FASTA_TRANSCRIPTS" ]; then
	echo "FASTA transcriptome file does not exist!" >&2
	exit 4
fi

# Check number of threads and set 1 as default value
if [ -z $THREADS ]; then
	THREADS=1
fi

# Check output
if [ -z "$OUTPUT" ]; then
	echo "Output directory must be specified!" >&2
	exit 5
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$OUTPUT")" ]; then
	echo "Output directory is not writable!" >&2
	exit 6
fi

#### Counting ####
SAMPLE_NAME=$(basename "$INPUT_BAM" ".bam")
SUFF="_sa.txt"
OUTPUT_NAME=$SAMPLE_NAME$SUFF

TEMP_DIR="TMP"

sudo docker run -v `pwd`:`pwd` -w `pwd` combinelab/salmon salmon quant -t "$FASTA_TRANSCRIPTS" -l A -a "$INPUT_BAM"  -p "$THREADS" -o "$OUTPUT"/"$TEMP_DIR"

OUTPUT_FILE="$OUTPUT"/"$TEMP_DIR"/quant.sf

# Check output file
if [ ! -f "$OUTPUT_FILE" ]; then
	echo "Unable to find output file!" >&2
	exit 7
fi

# Move output file from tmp directory to output directory
sudo chmod -R 777 "$OUTPUT"/"$TEMP_DIR"
mv "$OUTPUT_FILE" "$OUTPUT"/"$OUTPUT_NAME"

# Removing items of tmp directory
if [ -d "$OUTPUT"/"$TEMP_DIR" ]; then
	rm -rf "$OUTPUT"/"$TEMP_DIR"
fi
