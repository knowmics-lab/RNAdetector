#!/bin/bash

##############################################################################
# Options:
#   -n Package Name
# 	-u Package URL
#	-m Package MD5
##############################################################################
while getopts ":n:u:m:" opt; do
    case $opt in
    n) NAME=$OPTARG ;;
    u) URL=$OPTARG ;;
    m) MD5=$OPTARG ;;
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
if [ -z "$NAME" ]; then
    echo "Package Name is required."
    exit 4
fi

if [ -z "$URL" ]; then
    echo "Package URL is required."
    exit 3
fi

if [ -z "$MD5" ]; then
    echo "Package MD5 is required."
    exit 4
fi

REF_DIR="/rnadetector/ws/storage/app/references/"
FILENAME="${REF_DIR}/${NAME}.tar.bz2"
MD5_FILE="${REF_DIR}/${NAME}.tar.bz2.md5"

if [ -f "$FILENAME" ]; then
    echo "Package file already exists!"
else
    echo "Downloading package..."
    if ! curl -fSL "$URL" -o "$FILENAME"; then
        echo "Unable to download the package"
        exit 5
    fi
fi

if [ ! -f "$MD5_FILE" ]; then
    echo "Downloading package checksum file..."
    if ! curl -fSL "$MD5" -o "$MD5_FILE"; then
        echo "Unable to download MD5 checksum of the package"
        exit 6
    fi
fi

CURR_PWD=$(pwd)
cd $REF_DIR

echo "Checking package integrity..."
if ! md5sum -c "$MD5_FILE"; then
    echo "Checksum control failed."
    rm $FILENAME
    rm $MD5_FILE
    exit 7
fi

echo "Installing package..."
if ! php /rnadetector/ws/artisan reference:import "$NAME"; then
    echo "Unable to install package."
    exit 8
fi

rm $FILENAME
rm $MD5_FILE

echo "Package installed!"
cd "$CURR_PWD"
