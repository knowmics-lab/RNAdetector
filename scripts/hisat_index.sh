#!/bin/bash

##############################################################################
# Options:
# 	-f FASTA input genome file
#   -p STR	The basename of the index files to write
# 	-t NUMBER OF THREADS
##############################################################################
while getopts ":f:p:t:" opt; do
	case $opt in
	f) FASTA_FILE=$OPTARG ;;
	p) PREFIX_OUTPUT=$OPTARG ;;
	t) THREADS=$OPTARG ;;
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
# Check number of threads and set 1 as default value
if [ -z "$THREADS" ]; then
	THREADS=1
fi

# Check input files
if [ -z "$FASTA_FILE" ] || [ ! -f "$FASTA_FILE" ]; then
	echo "Input file does not exist!"
	exit 3
fi

# Check output
if [ -z "$PREFIX_OUTPUT" ]; then
	echo "Output prefix must be specified!"
	exit 4
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$PREFIX_OUTPUT")" ]; then
	echo "Output directory is not writable!"
	exit 5
fi

#### Genome indexing ####
if ! hisat2-build -p "$THREADS" "$FASTA_FILE" "$PREFIX_OUTPUT"; then
	echo "An error occurred during hisat2-build execution!"
	exit 6
fi

chmod -R 777 "$(dirname "$PREFIX_OUTPUT")"
