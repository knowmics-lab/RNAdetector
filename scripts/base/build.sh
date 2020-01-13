#!/usr/bin/env bash

if ! git clone https://github.com/alessandrolaferlita/RNAdetector.git; then
    exit 1
fi
tar -zcvf repo.tar.gz --owner=0 --group=0 ./RNAdetector
rm -rf ./RNAdetector
docker build -t alaimos/ubuntu-private:RNAdetector.v1.0 .
rm repo.tar.gz
