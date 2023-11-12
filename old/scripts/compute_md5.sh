#!/bin/bash

##############################################################################
# Options:
# 	-i input file
##############################################################################

exit_abnormal() {
  echo "$1" 1>&2
  # shellcheck disable=SC2086
  exit $2
}

while getopts ":i:" opt; do
  case $opt in
  i) INPUT_FILE=$OPTARG ;;
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
{ [ -z "$INPUT_FILE" ] || [ ! -f "$INPUT_FILE" ]; } && exit_abnormal "Input file does not exist!" 3

OUTPUT_FILE="${INPUT_FILE}.md5"

md5sum "$INPUT_FILE" > "$OUTPUT_FILE"  || exit_abnormal "Unable to compute MD5" 6
chmod 777 "$OUTPUT_FILE"
