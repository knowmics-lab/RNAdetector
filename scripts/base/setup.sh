#!/usr/bin/env bash

# Create Web Service Directory
mkdir /rnadetector/

# Create and fill temporary directory
mkdir /rnadetector/tmp/
cd /rnadetector/tmp/
tar -zxvf /repo.tar.gz
mv /rnadetector/tmp/RNAdetector/WS/ /rnadetector/ws/
mv /rnadetector/tmp/RNAdetector/scripts/ /rnadetector/scripts/
rm -rf /rnadetector/scripts/base/
rm /repo.tar.gz

# Install trim_galore
cd /rnadetector/tmp/
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.tar.gz -o trim_galore.tar.gz
tar -zxvf trim_galore.tar.gz
cp TrimGalore-0.6.5/trim_galore /usr/local/bin/
chmod 755 /usr/local/bin/

# Install the web service
cd /rnadetector/ws/
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
sed -i 's/post_max_size \= .M/post_max_size \= 100G/g'             /etc/php/*/fpm/php.ini
sed -i 's/upload_max_filesize \= .M/upload_max_filesize \= 100G/g' /etc/php/*/fpm/php.ini

# Redirect NGINX and PHP log to docker stdout and stderr
if [ -f "/var/log/nginx/access.log" ]; then
    rm /var/log/nginx/access.log
fi
ln -s /dev/stdout /var/log/nginx/access.log

if [ -f "/var/log/nginx/error.log" ]; then
    rm /var/log/nginx/error.log
fi
ln -s /dev/stderr /var/log/nginx/error.log 

if [ -f "/var/log/php7.2-fpm.log" ]; then
    rm /var/log/php7.2-fpm.log
fi
ln -s /dev/stderr /var/log/php7.2-fpm.log

# Set folder permission
chmod 755 /usr/local/bin/bootstrap.sh
chmod 755 /usr/local/bin/CIRI.pl
chmod 755 /rnadetector/scripts/*
chmod -R 777 /rnadetector/ws/bootstrap/cache
chmod -R 777 /rnadetector/ws/storage
chmod 755 /genkey.sh
chmod 755 /import_reference.sh
