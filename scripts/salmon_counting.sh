#!/bin/bash

##############################################################################
# Options:
#   -i Indexed trascriptome folder (input)
#   -f FIRST INPUT FASTQ (trimmed FASTQ file)
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED) (trimmed FASTQ file)
# 	-t NUMBER OF THREADS
#		-o OUTPUT file
#   -x map file
#   -h HARMONIZED gene counts
#   -n HARMONIZED transcripts counts
##############################################################################
while getopts ":i:f:s:t:o:h:n:x:" opt; do
    case $opt in
    i) INDEXED_FASTA=$OPTARG ;;
    f) INPUT_1=$OPTARG ;;
    s) INPUT_2=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    h) HARMONIZED=$OPTARG ;;
    n) HARMONIZED_TX=$OPTARG ;;
    x) MAP_FILE=$OPTARG ;;
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

# Check output
if [ -z "$OUTPUT" ]; then
    echo "Output file must be specified!"
    exit 5
fi

# Check if output directory is writable
if [ ! -w "$(dirname "$OUTPUT")" ]; then
    echo "Output directory is not writable!"
    exit 6
fi

# Check indexed fasta file
if [ ! -d "$INDEXED_FASTA" ]; then
    echo "Indexed trascriptome does not exist!"
    exit 7
fi

#### Counting ####
OUTPUT_DIR=$(dirname "$OUTPUT")
TEMP_DIR="$OUTPUT_DIR/TMP"

if [ $PAIRED = "true" ]; then
    if ! salmon quant -i "$INDEXED_FASTA" -l A -1 "$INPUT_1" -2 "$INPUT_2" --validateMappings -p $THREADS -o "$TEMP_DIR"; then
        echo "An error occurred during salmon quant execution!"
        exit 9
    fi
else
    if ! salmon quant -i "$INDEXED_FASTA" -l A -r "$INPUT_1" --validateMappings -p $THREADS -o "$TEMP_DIR"; then
        echo "An error occurred during salmon quant execution!"
        exit 9
    fi
fi

OUTPUT_FILE="$TEMP_DIR/quant.sf"

# Check output file
if [ ! -f "$OUTPUT_FILE" ]; then
    echo "Unable to find output file!"
    exit 8
fi

# Move output file from tmp directory to output directory
mv "$OUTPUT_FILE" "$OUTPUT"
chmod 777 "$OUTPUT"

# Removing items of tmp directory
if [ -d "$TEMP_DIR" ]; then
    rm -rf "$TEMP_DIR"
fi

if [ ! -z "$HARMONIZED" ]; then
    if [ -z "$HARMONIZED_TX" ]; then
        echo "Harmonized transcripts output file is required"
        exit 10
    fi
    CURR_DIR=$(pwd)
    SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
    cd $CURR_DIR
    if [ ! -z "$MAP_FILE" ] && [ -f "$MAP_FILE" ]; then
        if ! Rscript "${SCRIPT_PATH}/harmonize.R" -i "$OUTPUT" -a "salmon" -o "$HARMONIZED" -t "$HARMONIZED_TX" -m "$MAP_FILE"; then
            echo "Unable to harmonize output file"
            exit 11
        fi
    else
        if ! Rscript "${SCRIPT_PATH}/harmonize.R" -i "$OUTPUT" -a "salmon" -o "$HARMONIZED" -t "$HARMONIZED_TX"; then
            echo "Unable to harmonize output file"
            exit 11
        fi
    fi
    chmod 777 "$HARMONIZED"
    chmod 777 "$HARMONIZED_TX"
fi
