#!/bin/bash

##############################################################################
# Options:
#   -r FASTA transcripts files
#   -i input BAM file
# 	-t NUMBER OF THREADS
#		-o OUTPUT file
##############################################################################
while getopts ":r:i:t:o:" opt; do
    case $opt in
        r) FASTA_TRANSCRIPTS=$OPTARG ;;
        i) INPUT_BAM=$OPTARG ;;
        t) THREADS=$OPTARG ;;
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
# Check input files
if [ -z "$INPUT_BAM" ] || [ ! -f "$INPUT_BAM" ]; then
    echo "Input file does not exist!"
    exit 3
fi

# Check check fasta transcriptome
if [ -z "$FASTA_TRANSCRIPTS" ] || [ ! -f "$FASTA_TRANSCRIPTS" ]; then
    echo "FASTA transcriptome file does not exist!"
    exit 4
fi

# Check number of threads and set 1 as default value
if [ -z "$THREADS" ]; then
    THREADS=1
fi

# Check output
if [ -z "$OUTPUT" ]; then
    echo "Output directory must be specified!"
    exit 5
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$OUTPUT")" ]; then
    echo "Output directory is not writable!"
    exit 6
fi

#### Counting ####
OUTPUT_DIR=$(dirname "$OUTPUT")
TEMP_DIR="$OUTPUT_DIR/TMP"

salmon quant -t "$FASTA_TRANSCRIPTS" -l A -a "$INPUT_BAM" -p "$THREADS" -o "$TEMP_DIR"

OUTPUT_FILE="$TEMP_DIR/quant.sf"

# Check output file
if [ ! -f "$OUTPUT_FILE" ]; then
    echo "Unable to find output file!"
    exit 7
fi

# Move output file from tmp directory to output directory
mv "$OUTPUT_FILE" "$OUTPUT"
chmod 777 "$OUTPUT"

# Removing items of tmp directory
if [ -d "$TEMP_DIR" ]; then
    rm -rf "$TEMP_DIR"
fi
