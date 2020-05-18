#!/bin/bash

INPUT="${1}"

if [ -z $INPUT ]; then
  echo "Input file is required"
  exit 1
fi

if [ ! -f $INPUT ]; then
  echo "Input file does not exist."
  exit 2
fi

MIME=$(file -i -0 "$INPUT" | cut -f 2 -d " " | cut -f 1 -d ";")
EXTENSION="${INPUT##*.}"
FILENAME="$(dirname "$INPUT")/$(basename "$INPUT" ".${EXTENSION}")"

COMPRESSED=false
if [ "$MIME" = "application/gzip" ]; then
  echo "Extracting input file"
  if ! gunzip -k -v -c "$INPUT" 2>/dev/null >"$FILENAME" then
    echo "Unable to extract input file."
    exit 5
  fi
  COMPRESSED=true
elif [ "$MIME" = "application/x-bzip2" ]; then
  echo "Extracting input file"
  if ! bunzip2 -k -v -c "$INPUT" 2>/dev/null >"$FILENAME" then
    echo "Unable to extract input file."
    exit 5
  fi
  COMPRESSED=true
fi

if [ "$COMPRESSED" = "true" ]; then
  if [ ! -f "$FILENAME" ]; then
    echo "Unable to find extracted file"
    exit 3
  fi
  exit 4
fi
