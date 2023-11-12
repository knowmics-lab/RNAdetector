#!/bin/bash

[[ $DEBUG == true ]] && set -x

/usr/bin/mysqld_safe > /dev/null 2>&1 &

RET=1
while [[ RET -ne 0 ]]; do
    echo "Waiting for MySQL service startup..."
    sleep 5
    mysql -uroot -e "status" > /dev/null 2>&1
    RET=$?
done

php /rnadetector/ws/artisan migrate:fresh --seed --force || exit 100
php /rnadetector/ws/artisan first:boot || exit 101

mysqladmin -uroot shutdown
