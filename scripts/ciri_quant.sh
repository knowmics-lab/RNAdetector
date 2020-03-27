#!/bin/bash

##############################################################################
# Options:
#   -t number of threads (default: 1)
# 	-f FIRST INPUT FASTQ (trimmed FASTQ file)
# 	-s SECOND INPUT FASTQ (FOR PAIRED) (trimmed FASTQ file)
#   -c YAML-formated config file
#   -b Back-Spliced Junction Site in BED format
#   -o output file (required)
#   -x map file
#   -h HARMONIZED output file
##############################################################################
while getopts ":t:f:s:c:b:o:h:" opt; do
	case $opt in
	t) THREADS=$OPTARG ;;
	f) INPUT_1=$OPTARG ;;
	s) INPUT_2=$OPTARG ;;
	c) CONFIG_FILE=$OPTARG ;;
	b) BED_FILE=$OPTARG ;;
	o) OUTPUT=$OPTARG ;;
	x) MAP_FILE=$OPTARG ;;
	h) HARMONIZED=$OPTARG ;;
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
	echo "First input file does not exist!"
	exit 3
fi

if [ -z "$INPUT_2" ] || [ ! -f "$INPUT_2" ]; then
	echo "Second input file does not exist!"
	exit 4
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
	echo "Output directory must be specified!"
	exit 7
fi

OUT_DIR="$(dirname "$OUTPUT")"
OUT_PREFIX="$(basename "$OUTPUT" ".gtf")"

# Check if output directory is writable
if [ ! -w "$OUT_DIR" ]; then
	echo "Output directory is not writable!"
	exit 8
fi

#### Run CIRIquant ####
if ! CIRIquant -t "$THREADS" --bed "$BED_FILE" -o "$OUT_DIR" -p "$OUT_PREFIX" --config "$CONFIG_FILE" -1 "$INPUT_1" -2 "$INPUT_2"; then
	echo "An error occurred during CIRIquant execution!"
	exit 9
fi

# Check SAM file
if [ ! -f "$OUTPUT" ]; then
	echo "Unable to find CIRIquant output file!"
	exit 10
fi

chmod 777 "$OUTPUT"

if [ ! -z "$HARMONIZED" ]; then
	CURR_DIR=$(pwd)
	SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
	cd $CURR_DIR
	if [ ! -z "$MAP_FILE" ] && [ -f "$MAP_FILE" ]; then
		if ! Rscript "${SCRIPT_PATH}/harmonize.R" -i "$OUTPUT" -a "ciriquant" -o "$HARMONIZED" -m "$MAP_FILE"; then
			echo "Unable to harmonize output file"
			exit 11
		fi
	else
		if ! Rscript "${SCRIPT_PATH}/harmonize.R" -i "$OUTPUT" -a "ciriquant" -o "$HARMONIZED"; then
			echo "Unable to harmonize output file"
			exit 11
		fi
	fi
	chmod 777 "$HARMONIZED"
fi
