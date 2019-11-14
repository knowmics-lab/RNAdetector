#!/bin/bash

##############################################################################
# Options:
# 	-f FASTA input genome file
#   -p STR	The basename of the index files to write
##############################################################################
while getopts ":f:p:" opt; do
	case $opt in
		f ) FASTA_FILE=$OPTARG ;;
    p ) PREFIX_OUTPUT=$OPTARG;;
		\?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
		: ) echo "Option -$OPTARG requires an argument." >&2; exit 2;;
	esac
done

#### Check parameters ####

# Check input files
if [ -z "$FASTA_FILE" ] || [ ! -f "$FASTA_FILE" ]; then
	echo "Input file does not exist!" >&2
	exit 3
fi

# Check output
if [ -z  "$PREFIX_OUTPUT" ]; then
	echo "Output prefix must be specified!" >&2
	exit 4
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$PREFIX_OUTPUT")" ]; then
	echo "Output directory is not writable!" >&2
	exit 5
fi

#### Genome indexing ####
bowtie2-build  "$FASTA_FILE" "$PREFIX_OUTPUT"
