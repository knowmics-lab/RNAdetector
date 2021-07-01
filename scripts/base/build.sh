#!/usr/bin/env bash

if ! git clone https://github.com/alessandrolaferlita/RNAdetector.git; then
    exit 1
fi
if [ ! -z "$1" ]; then
    cd ./RNAdetector
    git checkout "$1"
    cd ..
fi
tar -zcvf repo.tar.gz --owner=0 --group=0 ./RNAdetector
rm -rf ./RNAdetector
docker build -t alaimos/rnadetector:v0.0.3 .
rm repo.tar.gz
