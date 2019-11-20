#!/usr/bin/env bash

git clone https://github.com/alessandrolaferlita/RNAdetector.git
tar -zcvf repo.tar.gz ./RNAdetector
rm -rf ./RNAdetector
docker build -t alaimos/ubuntu-base .
rm repo.tar.gz
