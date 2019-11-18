#!/bin/bash

##############################################################################
# Options:
# 	-r FASTA transcripts files
#   -i Indexed trascriptome (output folder)
##############################################################################
while getopts ":r:i:" opt; do
	case $opt in
		r ) FASTA=$OPTARG ;;
		i ) INDEXED_FASTA=$OPTARG ;;
		\?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
		: ) echo "Option -$OPTARG requires an argument." >&2; exit 2;;
	esac
done

#### Check parameters ####
# Check FASTA file
if [ -z "$FASTA" ] || [ ! -f "$FASTA" ]; then
	echo "FASTA file with transcripts does not exist!"
	exit 3
fi

#### Indexed transcriptome ####
salmon index -t "$FASTA" -i "$INDEXED_FASTA" -k 31

if [ ! -d "$INDEXED_FASTA" ]; then
	echo "Indexed trascriptome does not exist!"
	exit 4
fi
