#!/bin/bash

[[ $DEBUG == true ]] && set -x

DB_NAME="rnadetector"
DB_USER="rnadetector"
DB_PASS="secret"
MYSQL_CHARSET="utf8"
MYSQL_COLLATION="utf8_unicode_ci"

/usr/bin/mysqld_safe > /dev/null 2>&1 &

RET=1
while [[ RET -ne 0 ]]; do
    echo "Waiting for MySQL service startup..."
    sleep 5
    mysql -uroot -e "status" > /dev/null 2>&1
    RET=$?
done

[[ $DEBUG == true ]] && mysql -uroot -e "CREATE USER 'admin'@'%' IDENTIFIED BY '${DB_PASS}'"
[[ $DEBUG == true ]] && mysql -uroot -e "GRANT ALL PRIVILEGES ON *.* TO 'admin'@'%' WITH GRANT OPTION"

mysql -uroot -e "CREATE USER '${DB_USER}'@'%' IDENTIFIED BY '${DB_PASS}'"
mysql -uroot -e "GRANT USAGE ON *.* TO '${DB_USER}'@'%'"
mysql -uroot -e "CREATE DATABASE IF NOT EXISTS \`${DB_NAME}\` DEFAULT CHARACTER SET \`${MYSQL_CHARSET}\` COLLATE \`${MYSQL_COLLATION}\`"
mysql -uroot -e "GRANT ALL PRIVILEGES ON ${DB_NAME}.* TO '${DB_USER}'@'%'"

php /rnadetector/ws/artisan migrate --seed --force || exit 100
php /rnadetector/ws/artisan first:boot || exit 101

mysqladmin -uroot shutdown
