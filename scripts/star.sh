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

{ [ -z "$INPUT_1" ] || [ ! -f "$INPUT_1" ]; } && exit_abnormal "Input file does not exist!" 3

if [ -z "$INPUT_2" ]; then
  PAIRED=false
elif [ ! -f "$INPUT_2" ]; then
  exit_abnormal "Second input file does not exist!" 4
else
  PAIRED=true
fi

# Check number of threads and set 1 as default value
[ -z "$THREADS" ] && THREADS=1

# Check output
[ -z "$OUTPUT" ] && exit_abnormal "Output file must be specified!" 5
[ ! -w "$(dirname "$OUTPUT")" ] && exit_abnormal "Output directory is not writable!" 6

REFERENCE_DIR="${REF_GENOME}_star"

#### Alignment ####
echo "Detecting maximum read size"
MAX_SIZE=$(awk 'BEGIN {max=0} NR%4 == 2 {if (length($0)>max) {max=length($0)}} END {print max - 1}' "$INPUT_1")
echo "Setting sjdbOverhang to ${MAX_SIZE}."

TEMP_DIR="$(dirname "$OUTPUT")/star_tmp/"

([ ! -d "$TEMP_DIR" ] && mkdir -p "$TEMP_DIR") || exit_abnormal "Unable to create temp directory" 7

if [ $PAIRED = "true" ]; then
  STAR --runThreadN "$THREADS" --runMode alignReads \
    --genomeDir "$REFERENCE_DIR" --sjdbGTFfile "$GTF_FILE" \
    --sjdbOverhang "$MAX_SIZE" --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
    --readFilesIn "$INPUT_1" "$INPUT_2" \
    --outFileNamePrefix "$TEMP_DIR" || exit_abnormal "An error occurred during STAR execution!" 8
else
  STAR --runThreadN "$THREADS" --runMode alignReads \
    --genomeDir "$REFERENCE_DIR" --sjdbGTFfile "$GTF_FILE" \
    --sjdbOverhang "$MAX_SIZE" --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
    --readFilesIn "$INPUT_1" \
    --outFileNamePrefix "$TEMP_DIR" || exit_abnormal "An error occurred during STAR execution!" 8
fi

[ ! -f "$TEMP_DIR/Aligned.sortedByCoord.out.bam" ] && exit_abnormal "Unable to find STAR output file!" 9

[ -f "$TEMP_DIR/Log.final.out" ] && cat "$TEMP_DIR/Log.final.out"

mv "$TEMP_DIR/Aligned.sortedByCoord.out.bam" "$OUTPUT" || exit_abnormal "Unable to move output file!" 10

[ ! -f "$OUTPUT" ] && exit_abnormal "Unable to find output file!" 10

rm -r "$TEMP_DIR"

samtools index "$OUTPUT" || exit_abnormal "Unable to write index file!" 11

echo "Computing BAM coverage"
bamCoverage -b "$OUTPUT" -o "$OUTPUT.coverage.bw" || exit_abnormal "Unable to compute coverate!" 12

chmod 777 "$OUTPUT"
chmod 777 "$OUTPUT.bai"
chmod 777 "$OUTPUT.coverage.bw"
