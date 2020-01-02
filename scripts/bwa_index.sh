#!/bin/bash

##############################################################################
# Options:
# 	-f FASTA input genome file
#   -p STR	Prefix of the output database [same as db filename]
##############################################################################
while getopts ":f:p:a:" opt; do
	case $opt in
	f) FASTA_FILE=$OPTARG ;;
	p) PREFIX_OUTPUT=$OPTARG ;;
	a) ALGORITHM=$OPTARG ;; # -a option is now ignored
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
if ! bwa index -p "$PREFIX_OUTPUT" -a bwtsw "$FASTA_FILE"; then
	echo "An error occurred during bwa index execution!"
	exit 6
fi

chmod -R 777 "$(dirname "$PREFIX_OUTPUT")"
