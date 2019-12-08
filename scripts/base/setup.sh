#!/usr/bin/env bash

# Create Web Service Directory
mkdir /rnadetector/

# Create and fill temporary directory
mkdir /rnadetector/tmp/
cd /rnadetector/tmp/ || exit 100
tar -zxvf /repo.tar.gz
mv /rnadetector/tmp/RNAdetector/WS/ /rnadetector/ws/
mv /rnadetector/tmp/RNAdetector/scripts/ /rnadetector/scripts/
rm -rf /rnadetector/scripts/base/
rm /repo.tar.gz

# Install trim_galore
cd /rnadetector/tmp/ || exit 100
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.tar.gz -o trim_galore.tar.gz
tar -zxvf trim_galore.tar.gz
cp TrimGalore-0.6.5/trim_galore /usr/local/bin/
chmod 755 /usr/local/bin/

# Install latest version of salmon
cd /rnadetector/tmp || exit 100
curl -fsSL https://github.com/COMBINE-lab/salmon/releases/download/v1.0.0/salmon-1.0.0_linux_x86_64.tar.gz -o salmon.tar.gz
tar -zxvf salmon.tar.gz
mv salmon-latest_linux_x86_64/ /opt/salmon/
ln -s /opt/salmon/bin/salmon /usr/bin/salmon

# Install bbmap for re-pair utility
curl -fsSL https://sourceforge.net/projects/bbmap/files/latest/download -o bbmap.tar.gz
tar -zxvf bbmap.tar.gz --directory=/opt/
chmod 755 /opt/bbmap/*
rm bbmap.tar.gz

# Install stringtie
curl -fsSL http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.0.5.Linux_x86_64.tar.gz -o stringtie.tar.gz
tar -xzvf stringtie.tar.gz
cp stringtie-2.0.5.Linux_x86_64/stringtie /usr/local/bin/stringtie
cp stringtie-2.0.5.Linux_x86_64/prepDE.py /usr/local/bin/prepDE.py
rm -rf stringtie-2.0.5.Linux_x86_64/
rm stringtie.tar.gz

# Install latest version of htseq-count
pip install HTSeq

# Install latest version of fastq-pair
cd /rnadetector/tmp || exit 100
git clone https://github.com/linsalrob/fastq-pair.git
cd fastq-pair/ || exit 100
mkdir build && cd build && cmake .. && make && make install

# Install the web service
cd /rnadetector/ws/ || exit 100
mv .env.docker .env
composer install --optimize-autoloader --no-dev
php artisan key:generate
php artisan storage:link

# Remove temporary directory
rm -rf /rnadetector/tmp

# Copy configuration files
rm /etc/nginx/sites-available/default
mv /nginx.conf /etc/nginx/sites-available/default
mv /worker.conf /etc/supervisor/conf.d/worker.conf
sed -i 's/post_max_size \= .M/post_max_size \= 100G/g' /etc/php/*/fpm/php.ini
sed -i 's/upload_max_filesize \= .M/upload_max_filesize \= 100G/g' /etc/php/*/fpm/php.ini

# Redirect NGINX and PHP log to docker stdout and stderr
if [ -f "/var/log/nginx/access.log" ]; then
    rm /var/log/nginx/access.log
fi
ln -s /dev/stdout /var/log/nginx/access.log

if [ -f "/var/log/nginx/error.log" ]; then
    rm /var/log/nginx/error.log
fi
ln -s /dev/stdout /var/log/nginx/error.log

if [ -f "/var/log/php7.3-fpm.log" ]; then
    rm /var/log/php7.3-fpm.log
fi
ln -s /dev/stdout /var/log/php7.3-fpm.log

# Set folder permission
chmod 755 /usr/local/bin/bootstrap.sh
chmod 755 /usr/local/bin/CIRI1.pl
chmod 755 /usr/local/bin/CIRI2.pl
chmod 755 /rnadetector/scripts/*
chmod -R 777 /rnadetector/ws/bootstrap/cache
chmod -R 777 /rnadetector/ws/storage
chmod 755 /genkey.sh
chmod 755 /import_reference.sh
