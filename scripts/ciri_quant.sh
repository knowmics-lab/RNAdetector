#!/bin/bash

##############################################################################
# Options:
#   -t number of threads (default: 1)
# 	-f FIRST INPUT FASTQ (trimmed FASTQ file)
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED) (trimmed FASTQ file)
#   -c YAML-formated config file
#   -b Back-Spliced Junction Site in BED format
#   -o output (required)
##############################################################################
while getopts ":t:f:s:c:b:o:" opt; do
	case $opt in
	t) THREADS=$OPTARG ;;
	f) INPUT_1=$OPTARG ;;
	s) INPUT_2=$OPTARG ;;
	c) CONFIG_FILE=$OPTARG ;;
	b) BED_FILE=$OPTARG ;;
	o) OUTPUT=$OPTARG ;;
	\?)
		echo "Invalid option: -$OPTARG"
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument."
		exit 2
		;;
	esac
done

#### Check parameters ####
#Check input files
if [ -z "$INPUT_1" ] || [ ! -f "$INPUT_1" ]; then
	echo "Input file does not exist!"
	exit 3
fi

# Control sequencing strategy "single end" o "paired end"
if [ -z "$INPUT_2" ]; then
	PAIRED=false
elif [ ! -f "$INPUT_2" ]; then
	echo "Second input file does not exist!"
	exit 4
else
	PAIRED=true
fi

# Check number of threads and set 1 as default value
if [ -z "$THREADS" ]; then
	THREADS=1
fi

# Check config file
if [ -z "$CONFIG_FILE" ] || [ ! -f "$CONFIG_FILE" ]; then
	echo "Configuration file does not exist!"
	exit 5
fi

# Check BED file
if [ -z "$BED_FILE" ] || [ ! -f "$BED_FILE" ]; then
	echo "BED file with Back-Spliced Junction Site does not exist!"
	exit 6
fi

# Check output
if [ -z "$OUTPUT" ]; then
	echo "Output file must be specified!"
	exit 7
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$OUTPUT")" ]; then
	echo "Output directory is not writable!"
	exit 8
fi

#### Run CIRIquant ####
if ! CIRIquant -t "$THREADS" --bed "$BED_FILE" -o "$OUTPUT" --config "$CONFIG_FILE" -1 "$INPUT_1" -2 "$INPUT_2"; then
	echo "An error occurred during CIRIquant execution!"
	exit 9
fi

# Check output file
if [ ! -f "$OUTPUT" ]; then
	echo "Unable to find output file!"
	exit 10
fi
