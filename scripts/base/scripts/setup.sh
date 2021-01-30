#!/usr/bin/env bash

# Update the distro
DEBIAN_FRONTEND=noninteractive /usr/bin/apt update
DEBIAN_FRONTEND=noninteractive /usr/bin/apt dist-upgrade -y

# Create Web Service Directory
[ -d /rnadetector ] && mkdir /rnadetector/

# Create and fill temporary directory
[ -d /rnadetector/tmp/ ] && mkdir /rnadetector/tmp/
cd /rnadetector/tmp/ || exit 100
tar -zxvf /repo.tar.gz
mv /rnadetector/tmp/RNAdetector/WS/ /rnadetector/ws/
mv /rnadetector/tmp/RNAdetector/scripts/ /rnadetector/scripts/
rm -rf /rnadetector/scripts/base/
rm /repo.tar.gz

# Install the web service
cd /rnadetector/ws/ || exit 100
mv .env.docker .env
composer install --optimize-autoloader --no-dev
php artisan key:generate
php artisan storage:link
php artisan make:links

# Download MITHrIL index
ORGANISMS="hsa rno mmu cel"
for O in $ORGANISMS; do
    java -jar /rnadetector/scripts/resources/pathways/MITHrIL2.jar index -enrichment-evidence-type STRONG -organism "$O" -verbose
    java -jar /rnadetector/scripts/resources/pathways/MITHrIL2.jar index -enrichment-evidence-type WEAK -organism "$O" -verbose
    java -jar /rnadetector/scripts/resources/pathways/MITHrIL2.jar index -enrichment-evidence-type PREDICTION -organism "$O" -verbose
done

# Remove temporary directory
rm -rf /rnadetector/tmp

# Apply PHP configuration fixes
sed -i 's/post_max_size \= .M/post_max_size \= 200G/g' /etc/php/*/apache2/php.ini
sed -i 's/upload_max_filesize \= .M/upload_max_filesize \= 200G/g' /etc/php/*/apache2/php.ini
sed -i "s/;date.timezone =/date.timezone = Europe\/London/g" /etc/php/*/apache2/php.ini
sed -i "s/;date.timezone =/date.timezone = Europe\/London/g" /etc/php/*/cli/php.ini
sed -i "s/export APACHE_RUN_GROUP=www-data/export APACHE_RUN_GROUP=staff/" /etc/apache2/envvars


apply_configuration_fixes() {
  sed -i 's/^log_error/# log_error/' /etc/mysql/mysql.conf.d/mysqld.cnf
  sed -i 's/.*datadir.*/datadir = \/rnadetector\/ws\/storage\/app\/database/' /etc/mysql/mysql.conf.d/mysqld.cnf
  sed -i "s/.*bind-address.*/bind-address = 0.0.0.0/" /etc/mysql/my.cnf
  sed -i "s/.*bind-address.*/bind-address = 0.0.0.0/" /etc/mysql/mysql.conf.d/mysqld.cnf
  sed -i "s/user.*/user = www-data/" /etc/mysql/mysql.conf.d/mysqld.cnf
  cat >/etc/mysql/conf.d/mysql-skip-name-resolv.cnf <<EOF
[mysqld]
skip_name_resolve
EOF
}

remove_debian_system_maint_password() {
    sed 's/password = .*/password = /g' -i /etc/mysql/debian.cnf
}

apply_configuration_fixes
remove_debian_system_maint_password

# Set folder permission
chmod 755 /usr/local/bin/bootstrap.sh
chmod 755 /usr/local/bin/CIRI1.pl
chmod 755 /usr/local/bin/CIRI2.pl
chmod 755 /rnadetector/scripts/*
chmod -R 777 /rnadetector/ws/bootstrap/cache
chmod -R 777 /rnadetector/ws/storage
chmod 755 /genkey.sh
chmod 755 /import_reference.sh
