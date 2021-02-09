#!/bin/bash

##############################################################################
# Options:
# 	-f BAM input file
#   -s Enable sorting
#   -u Unsorted BAM file (required when -s is set)
#   -t the number of threads for sorting
##############################################################################

exit_abnormal() {
  echo "$1" 1>&2
  # shellcheck disable=SC2086
  exit $2
}

SORT=false
while getopts "sf:u:t:" opt; do
  case $opt in
  f) BAM=$OPTARG ;;
  s) SORT=true ;;
  u) UNSORTED=$OPTARG ;;
  t) THREADS=$OPTARG ;;
  \?)
    exit_abnormal "Invalid option: -$OPTARG" 1
    ;;
  :)
    exit_abnormal "Option -$OPTARG requires an argument." 2
    ;;
  esac
done

if [ "$SORT" = "true" ]; then
  { [ -z "$UNSORTED" ] || [ ! -f "$UNSORTED" ]; } && exit_abnormal "Unsorted BAM file does not exist!" 101
  [ -z "$THREADS" ] && THREADS=1
  echo "Sorting..."
  samtools sort --threads "$THREADS" "$UNSORTED" -o "$BAM" || exit_abnormal "Unable to sort BAM file!" 102
fi

{ [ -z "$BAM" ] || [ ! -f "$BAM" ]; } && exit_abnormal "BAM file does not exist!" 103

echo "Indexing BAM file..."
samtools index "$BAM" || exit_abnormal "Unable to write index file!" 104

echo "Computing BAM coverage..."
bamCoverage -b "$BAM" -o "$BAM.coverage.bw" || exit_abnormal "Unable to compute coverage!" 105

chmod 777 "$BAM"
chmod 777 "$BAM.bai"
chmod 777 "$BAM.coverage.bw"
