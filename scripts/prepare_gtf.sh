#!/bin/bash

##############################################################################
# Options:
# 	-f GTF input file
##############################################################################

exit_abnormal() {
  echo "$1" 1>&2
  # shellcheck disable=SC2086
  exit $2
}

while getopts ":f:" opt; do
  case $opt in
  f) GTF_FILE=$OPTARG ;;
  \?)
    exit_abnormal "Invalid option: -$OPTARG" 1
    ;;
  :)
    exit_abnormal "Option -$OPTARG requires an argument." 2
    ;;
  esac
done

#### Check parameters ####
# Check input files
{ [ -z "$GTF_FILE" ] || [ ! -f "$GTF_FILE" ]; } && exit_abnormal "Input file does not exist!" 3

OUTPUT_DIRECTORY=$(dirname "$GTF_FILE")
# Check if output directory is writable
[ ! -w "$OUTPUT_DIRECTORY" ] && exit_abnormal "Output directory is not writable!" 4

NAME=$(basename "$GTF_FILE" ".gtf")

gffread -E "$GTF_FILE" -o "$OUTPUT_DIRECTORY/$NAME.gff3" || exit_abnormal "Unable to convert GTF to GFF3" 5
bedtools sort -i "$OUTPUT_DIRECTORY/$NAME.gff3" | bgzip > "$OUTPUT_DIRECTORY/$NAME.gff3.gz"  || exit_abnormal "Unable to sort GFF3" 6
tabix -p gff "$OUTPUT_DIRECTORY/$NAME.gff3.gz" || exit_abnormal "Unable to index GFF3" 7
rm "$OUTPUT_DIRECTORY/$NAME.gff3" || exit_abnormal "Unable to remove temporary files" 8

chmod -R 777 "$OUTPUT_DIRECTORY/$NAME.gff3.gz"
chmod -R 777 "$OUTPUT_DIRECTORY/$NAME.gff3.gz.tbi"
