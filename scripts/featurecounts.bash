#!/bin/bash

##############################################################################
# Options:
# 	-a FILE GTF
#	-b INPUT BAM (BAM_PATH or ALIGNMENT_PATH)
# 	-t NUMBER OF THREADS
#	-o OUTPUT file
#   -x map file
#   -h HARMONIZED output file
##############################################################################
OTHER_ARGS=""
while getopts ":a:b:t:o:h:x:A:" opt; do
  case $opt in
  a) GTF_FILE=$OPTARG ;;
  b) INPUT_BAM=$OPTARG ;;
  t) THREADS=$OPTARG ;;
  o) OUTPUT=$OPTARG ;;
  h) HARMONIZED=$OPTARG ;;
  x) MAP_FILE=$OPTARG ;;
  A) OTHER_ARGS=$OPTARG ;;
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

exiterror() {
  echo "$1"
  exit $2
}

#### Check parameters ####
# Check GTF annotation files
if [ -z "$GTF_FILE" ] || [ ! -f "$GTF_FILE" ]; then
  echo "Annotation file does not exist!"
  exit 3
fi

# Check input files
if [ -z "$INPUT_BAM" ] || [ ! -f "$INPUT_BAM" ]; then
  echo "Input file does not exist!"
  exit 4
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

if [ -n "$OTHER_ARGS" ]; then
  echo "Processing with custom arguments: \"${OTHER_ARGS}\""
fi

#### Counting ####
# shellcheck disable=SC2086
if ! featureCounts -T $THREADS $OTHER_ARGS -a "$GTF_FILE" -o "$OUTPUT" "$INPUT_BAM"; then
  echo "An error occurred during featureCounts execution!"
  exit 8
fi

if [ ! -f "$OUTPUT" ]; then
  echo "Unable to find output file!"
  exit 7
fi

chmod 777 "$OUTPUT"

if [ -n "$HARMONIZED" ]; then
  CURR_DIR=$(pwd)
  SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
  cd "$CURR_DIR" || exit 9
  [ -n "$MAP_FILE" ] && [ -f "$MAP_FILE" ] && (Rscript "${SCRIPT_PATH}/harmonize.R" -i "$OUTPUT" -g "$GTF_FILE" -a "featurecounts" -o "$HARMONIZED" -m "$MAP_FILE" || exiterror "Unable to harmonize output file" 9)
  ([ -z "$MAP_FILE" ] || [ ! -f "$MAP_FILE" ]) && (Rscript "${SCRIPT_PATH}/harmonize.R" -i "$OUTPUT" -g "$GTF_FILE" -a "featurecounts" -o "$HARMONIZED" || exiterror "Unable to harmonize output file" 9)
  chmod 777 "$HARMONIZED"
fi
