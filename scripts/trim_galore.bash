#!/bin/bash

##############################################################################
# Options:
#   -q QUALITY (default 20)
#   -l LENGTH (default 14)
#   -f FIRST INPUT FASTQ
#   -s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED)
#   -t NUMBER OF THREADS
#   -o OUTPUT TRIMMED FASTQ
##############################################################################

HARD_TRIM=false
OTHER_ARGS=""
while getopts ":q:l:f:s:t:o:A:h" opt; do
  case $opt in
  q) QUALITY=$OPTARG ;;
  l) LENGTH=$OPTARG ;;
  f) INPUT_1=$OPTARG ;;
  s) INPUT_2=$OPTARG ;;
  t) THREADS=$OPTARG ;;
  o) OUTPUT=$OPTARG ;;
  h) HARD_TRIM=true ;;
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

#### Check parameters ####
# Quality control
if [ -z "$QUALITY" ]; then
  QUALITY=20
fi
# Control read LENGTH
if [ -z "$LENGTH" ]; then
  LENGTH=14
fi
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

if ((THREADS >= 8)); then
  echo "Setting number of threads to 7 since a number greater than 8 can reduce performance!"
  THREADS=7
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

#### Trimming and adaptors removing ####
INPUT_1_NAME=$(basename -- "$INPUT_1")
INPUT_1_EXTENSION="${INPUT_1_NAME##*.}"
INPUT_1_BASENAME=$(basename "$INPUT_1_NAME" ".$INPUT_1_EXTENSION")
if [ $PAIRED = "true" ]; then
  # shellcheck disable=SC2086
  if ! trim_galore -j $THREADS -q $QUALITY --paired -o "$OUTPUT" --dont_gzip --phred33 --length $LENGTH --no_report_file $OTHER_ARGS "$INPUT_1" "$INPUT_2"; then
    echo "An error occurred during trim_galore execution!"
    exit 10
  fi
  if [ $HARD_TRIM = "true" ]; then
    INPUT_2_NAME=$(basename -- "$INPUT_2")
    INPUT_2_EXTENSION="${INPUT_2_NAME##*.}"
    INPUT_2_BASENAME=$(basename "$INPUT_2_NAME" ".$INPUT_2_EXTENSION")
    # shellcheck disable=SC2086
    if ! trim_galore -j $THREADS --paired -o "$OUTPUT" --dont_gzip --hardtrim5 $LENGTH --no_report_file $OTHER_ARGS "$OUTPUT/${INPUT_1_BASENAME}_val_1.fq" "$OUTPUT/${INPUT_2_BASENAME}_val_2.fq"; then
      echo "Failed to run trim_galore"
      exit 7
    fi
    rm "$OUTPUT/${INPUT_1_BASENAME}_val_1.fq"
    rm "$OUTPUT/${INPUT_2_BASENAME}_val_2.fq"
    OUT_FILE=$(find "$OUTPUT" -name "${INPUT_1_BASENAME}_val_1.*bp_5prime.fq" | head -n 1)
    if [ ! -f "$OUT_FILE" ]; then
      echo "Unable to find first hard trimmed output file"
      exit 8
    fi
    mv "$OUT_FILE" "$OUTPUT/${INPUT_1_BASENAME}_val_1.fq"
    OUT_FILE=$(find "$OUTPUT" -name "${INPUT_2_BASENAME}_val_2.*bp_5prime.fq" | head -n 1)
    if [ ! -f "$OUT_FILE" ]; then
      echo "Unable to find second hard trimmed output file"
      exit 9
    fi
    mv "$OUT_FILE" "$OUTPUT/${INPUT_2_BASENAME}_val_2.fq"
  fi
else
  # shellcheck disable=SC2086
  if ! trim_galore -j $THREADS -q $QUALITY -o "$OUTPUT" --dont_gzip --phred33 --length $LENGTH --no_report_file $OTHER_ARGS "$INPUT_1"; then
    echo "An error occurred during trim_galore execution!"
    exit 10
  fi
  if [ $HARD_TRIM = "true" ]; then
    # shellcheck disable=SC2086
    if ! trim_galore -j $THREADS -o "$OUTPUT" --dont_gzip --hardtrim5 $LENGTH --no_report_file $OTHER_ARGS "$OUTPUT/${INPUT_1_BASENAME}_trimmed.fq"; then
      echo "Failed to run trim_galore"
      exit 7
    fi
    rm "$OUTPUT/${INPUT_1_BASENAME}_trimmed.fq"
    OUT_FILE=$(find "$OUTPUT" -name "${INPUT_1_BASENAME}_trimmed.*bp_5prime.fq" | head -n 1)
    if [ ! -f "$OUT_FILE" ]; then
      echo "Unable to find hard trimmed output file"
      exit 8
    fi
    mv "$OUT_FILE" "$OUTPUT/${INPUT_1_BASENAME}_trimmed.fq"
  fi
fi
