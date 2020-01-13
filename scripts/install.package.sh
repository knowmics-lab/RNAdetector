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

cleanup() {
    if [ -f "$FILENAME" ]; then
        rm "$FILENAME"
    fi
    if [ -f "$MD5_FILE" ]; then
        rm "$MD5_FILE"
    fi
}


if [ -f "$FILENAME" ]; then
    echo "Package file already exists!"
else
    echo "Downloading package..."
    if ! wget "$URL" -O "$FILENAME"; then
        echo "Unable to download the package"
        cleanup
        exit 5
    fi
fi

if [ ! -f "$MD5_FILE" ]; then
    echo "Downloading package checksum file..."
    if ! wget "$MD5" -O "$MD5_FILE"; then
        echo "Unable to download MD5 checksum of the package"
        cleanup
        exit 6
    fi
fi

CURR_PWD=$(pwd)
cd $REF_DIR

echo "Checking package integrity..."
if ! md5sum -c "$MD5_FILE"; then
    echo "Checksum control failed."
    cleanup
    exit 7
fi

echo "Installing package..."
if ! php /rnadetector/ws/artisan reference:import "$NAME"; then
    echo "Unable to install package."
    cleanup
    exit 8
fi

cleanup

echo "Package installed!"
cd "$CURR_PWD"
