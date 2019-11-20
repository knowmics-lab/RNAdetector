#!/usr/bin/env bash

mkdir /rnadetector/
mkdir /rnadetector/tmp/
cd /rnadetector/tmp/
tar -zxvf /repo.tar.gz
mv /rnadetector/tmp/RNAdetector/WS/ /rnadetector/ws/
mv /rnadetector/tmp/RNAdetector/scripts/ /rnadetector/scripts/
rm -rf /rnadetector/scripts/base/
rm /repo.tar.gz

# curl -o CIRIquant_v1.0.tar.gz "https://liquidtelecom.dl.sourceforge.net/project/ciri/CIRIquant/CIRIquant_v1.0.tar.gz"
# cd /rnadetector/tmp/CIRIquant
# python setup.py install

cd /rnadetector/tmp/
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.tar.gz -o trim_galore.tar.gz
tar -zxvf trim_galore.tar.gz
cp TrimGalore-0.6.5/trim_galore /usr/local/bin/
chmod 755 /usr/local/bin/

rm -rf /rnadetector/tmp

cd /rnadetector/ws/
mv .env.docker .env
composer install --no-dev
php artisan key:generate
mkdir -p /rnadetector/ws/storage/app/database/
touch /rnadetector/ws/storage/app/database/database.sqlite
mkdir -p /rnadetector/ws/storage/app/annotations/
mkdir -p /rnadetector/ws/storage/app/references/
php artisan migrate --seed
