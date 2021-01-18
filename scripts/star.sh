#!/bin/bash

##############################################################################
# Options:
# 	-a FILE GTF
# 	-g REFERENCE INDEXED GENOME FILE (basename of indexed genome)
# 	-t NUMBER OF THREADS
# 	-f FIRST INPUT FASTQ (trimmed FASTQ file)
# 	-s OPTIONAL SECOND INPUT FASTQ (FOR PAIRED) (trimmed FASTQ file)
# 	-o OUTPUT BAM FILE
##############################################################################

exit_abnormal() {
  echo "$1" 1>&2
  # shellcheck disable=SC2086
  exit $2
}

while getopts ":a:g:t:f:s:o:" opt; do
  case $opt in
  a) GTF_FILE=$OPTARG ;;
  g) REF_GENOME=$OPTARG ;;
  t) THREADS=$OPTARG ;;
  f) INPUT_1=$OPTARG ;;
  s) INPUT_2=$OPTARG ;;
  o) OUTPUT=$OPTARG ;;
  \?)
    exit_abnormal "Invalid option: -$OPTARG" 1
    ;;
  :)
    exit_abnormal "Option -$OPTARG requires an argument." 2
    ;;
  esac
done

#### Check parameters ####
#Check input files
{ [ -z "$GTF_FILE" ] || [ ! -f "$GTF_FILE" ]; } && exit_abnormal "Annotation file does not exist!" 3

{ [ -z "$INPUT_1" ] || [ ! -f "$INPUT_1" ]; } && exit_abnormal "Input file does not exist!" 4

if [ -z "$INPUT_2" ]; then
  PAIRED=false
elif [ ! -f "$INPUT_2" ]; then
  exit_abnormal "Second input file does not exist!" 5
else
  PAIRED=true
fi

# Check number of threads and set 1 as default value
[ -z "$THREADS" ] && THREADS=1

# Check output
[ -z "$OUTPUT" ] && exit_abnormal "Output file must be specified!" 6
[ ! -w "$(dirname "$OUTPUT")" ] && exit_abnormal "Output directory is not writable!" 7

REFERENCE_DIR="${REF_GENOME}_star"

#### Alignment ####
echo "Detecting maximum read size"
MAX_SIZE=$(awk 'BEGIN {max=0} NR%4 == 2 {if (length($0)>max) {max=length($0)}} END {print max - 1}' "$INPUT_1")
echo "Setting sjdbOverhang to ${MAX_SIZE}."

TEMP_DIR="$(dirname "$OUTPUT")/star_tmp/"

([ ! -d "$TEMP_DIR" ] && mkdir -p "$TEMP_DIR") || exit_abnormal "Unable to create temp directory" 8

if [ $PAIRED = "true" ]; then
  STAR --runThreadN "$THREADS" --runMode alignReads \
    --genomeDir "$REFERENCE_DIR" --sjdbGTFfile "$GTF_FILE" \
    --sjdbOverhang "$MAX_SIZE" --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
    --readFilesIn "$INPUT_1" "$INPUT_2" \
    --outFileNamePrefix "$TEMP_DIR" || exit_abnormal "An error occurred during STAR execution!" 9
else
  STAR --runThreadN "$THREADS" --runMode alignReads \
    --genomeDir "$REFERENCE_DIR" --sjdbGTFfile "$GTF_FILE" \
    --sjdbOverhang "$MAX_SIZE" --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
    --readFilesIn "$INPUT_1" \
    --outFileNamePrefix "$TEMP_DIR" || exit_abnormal "An error occurred during STAR execution!" 9
fi

[ ! -f "$TEMP_DIR/Aligned.sortedByCoord.out.bam" ] && exit_abnormal "Unable to find STAR output file!" 10

[ -f "$TEMP_DIR/Log.final.out" ] && cat "$TEMP_DIR/Log.final.out"

mv "$TEMP_DIR/Aligned.sortedByCoord.out.bam" "$OUTPUT" || exit_abnormal "Unable to move output file!" 11

[ ! -f "$OUTPUT" ] && exit_abnormal "Unable to find output file!" 12

rm -r "$TEMP_DIR"

bash /rnadetector/scripts/prepare_bam.sh -f "$OUTPUT" || exit_abnormal "Unable to prepare BAM file" "$?"
