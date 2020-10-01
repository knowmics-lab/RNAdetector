#!/usr/bin/env bash
set -e

[[ $DEBUG == true ]] && set -x

MYSQL_DATA_DIR="/rnadetector/ws/storage/app/database/"
MYSQL_USER="mysql"
MYSQL_RUN_DIR="/var/run/mysqld"
DB_NAME="rnadetector"
DB_USER="rnadetector"
DB_PASS="secret"
MYSQL_CHARSET="utf8"
MYSQL_COLLATION="utf8_unicode_ci"
DB_CREATED="false"

create_data_dir() {
    if [ ! -d ${MYSQL_DATA_DIR} ]; then
        mkdir -p ${MYSQL_DATA_DIR}
    fi
    chmod -R 0777 ${MYSQL_DATA_DIR}
    chown -R ${MYSQL_USER}:${MYSQL_USER} ${MYSQL_DATA_DIR}
}

create_run_dir() {
    if [ ! -d ${MYSQL_RUN_DIR} ]; then
        mkdir -p ${MYSQL_RUN_DIR}
    fi
    chmod -R 0755 ${MYSQL_RUN_DIR}
    chown -R ${MYSQL_USER}:root ${MYSQL_RUN_DIR}
    rm -rf ${MYSQL_RUN_DIR}/mysqld.sock.lock
}

initialize_mysql_database() {
    # initialize MySQL data directory
    if [ ! -d ${MYSQL_DATA_DIR}/mysql ]; then
        echo "Installing database..."
        mysqld --initialize-insecure --user=mysql >/dev/stdout 2>&1
        echo "Starting MySQL server..."
        /usr/bin/mysqld_safe >/dev/stdout 2>&1 &
        timeout=30
        echo -n "Waiting for database server to accept connections"
        while ! /usr/bin/mysqladmin -u root status >/dev/null 2>&1; do
            timeout=$(($timeout - 1))
            if [ $timeout -eq 0 ]; then
                echo -e "\nCould not connect to database server. Aborting..."
                exit 1
            fi
            echo -n "."
            sleep 1
        done
        echo
        echo "Creating debian-sys-maint user..."
        mysql -uroot -e "CREATE USER 'debian-sys-maint'@'localhost' IDENTIFIED BY '';"
        mysql -uroot -e "GRANT ALL PRIVILEGES on *.* TO 'debian-sys-maint'@'localhost' IDENTIFIED BY '' WITH GRANT OPTION;"
        /usr/bin/mysqladmin --defaults-file=/etc/mysql/debian.cnf shutdown
    fi
}

create_users_and_databases() {
    if [ ! -d ${MYSQL_DATA_DIR}/${DB_NAME} ]; then
        /usr/bin/mysqld_safe >/dev/stdout 2>&1 &
        timeout=30
        while ! /usr/bin/mysqladmin -u root status >/dev/null 2>&1; do
            timeout=$(($timeout - 1))
            if [ $timeout -eq 0 ]; then
                echo "Could not connect to mysql server. Aborting..."
                exit 1
            fi
            sleep 1
        done
        echo "Creating database \"${DB_NAME}\"..."
        mysql --defaults-file=/etc/mysql/debian.cnf \
            -e "CREATE DATABASE IF NOT EXISTS \`${DB_NAME}\` DEFAULT CHARACTER SET \`${MYSQL_CHARSET}\` COLLATE \`${MYSQL_COLLATION}\`;"
        echo "Granting access to database \"${DB_NAME}\" for user \"${DB_USER}\"..."
        mysql --defaults-file=/etc/mysql/debian.cnf \
            -e "GRANT ALL PRIVILEGES ON \`${DB_NAME}\`.* TO '${DB_USER}' IDENTIFIED BY '${DB_PASS}';"
        if php /rnadetector/ws/artisan migrate --seed --force; then
            export DB_CREATED="true"
        fi
        /usr/bin/mysqladmin --defaults-file=/etc/mysql/debian.cnf shutdown
    fi
}

if [ ! -d "/rnadetector/ws/storage/app/public/" ]; then
    mkdir -p "/rnadetector/ws/storage/app/public/"
fi
if [ ! -d "/rnadetector/ws/storage/app/annotations/" ]; then
    mkdir -p "/rnadetector/ws/storage/app/annotations/"
fi
if [ ! -d "/rnadetector/ws/storage/app/references/" ]; then
    mkdir -p "/rnadetector/ws/storage/app/references/"
fi
if [ ! -d "/rnadetector/ws/storage/app/tus_cache/" ]; then
    mkdir -p "/rnadetector/ws/storage/app/tus_cache/"
fi
if [ ! -d "/rnadetector/ws/storage/app/logs/" ]; then
    mkdir -p "/rnadetector/ws/storage/app/logs/"
fi
create_data_dir
create_run_dir
initialize_mysql_database
create_users_and_databases

chmod -R 777 "/rnadetector/ws/storage/"

/etc/init.d/nginx start
/etc/init.d/php7.3-fpm start
/etc/init.d/mysql start
/etc/init.d/supervisor start
if [ "$DB_CREATED" = "true" ]; then
    touch "${MYSQL_DATA_DIR}/ready"
fi

exec "$@"
