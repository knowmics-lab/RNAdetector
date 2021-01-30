#!/usr/bin/env bash
set -e

[[ $DEBUG == true ]] && set -x

MYSQL_DATA_DIR="/rnadetector/ws/storage/app/database/"
MYSQL_USER="www-data"
MYSQL_GROUP="staff"
MYSQL_RUN_DIR="/var/run/mysqld"
DB_NAME="rnadetector"

create_data_dir() {
  if [ ! -d ${MYSQL_DATA_DIR} ]; then
    mkdir -p ${MYSQL_DATA_DIR}
  fi
  chmod -R 0777 ${MYSQL_DATA_DIR}
  chown -R ${MYSQL_USER}:${MYSQL_GROUP} ${MYSQL_DATA_DIR}
}

create_run_dir() {
  if [ ! -d ${MYSQL_RUN_DIR} ]; then
    mkdir -p ${MYSQL_RUN_DIR}
  fi
  chmod -R 0775 ${MYSQL_RUN_DIR}
  chown -R ${MYSQL_USER}:${MYSQL_GROUP} ${MYSQL_RUN_DIR}
  if [ -e ${MYSQL_RUN_DIR}/mysqld.sock ]; then
    rm ${MYSQL_RUN_DIR}/mysqld.sock
  fi
  rm -rf ${MYSQL_RUN_DIR}/mysqld.sock.lock
}

initialize_mysql_database() {
  # initialize MySQL data directory
  if [ ! -d ${MYSQL_DATA_DIR}/mysql ]; then
    echo "Installing database..."
    mysqld --initialize-insecure
    if [ $? -ne 0 ]; then
      mysql_install_db >/dev/null 2>&1
    fi
  fi
  if [ ! -d ${MYSQL_DATA_DIR}/${DB_NAME} ]; then
    echo "Creating users..."
    if /usr/local/bin/create_mysql_users.sh; then
      export DB_CREATED="true"
    fi
  fi
}

[ ! -d "/rnadetector/ws/storage/app/public/" ] && mkdir -p "/rnadetector/ws/storage/app/public/"
[ ! -d "/rnadetector/ws/storage/app/annotations/" ] && mkdir -p "/rnadetector/ws/storage/app/annotations/"
[ ! -d "/rnadetector/ws/storage/app/references/" ] && mkdir -p "/rnadetector/ws/storage/app/references/"
[ ! -d "/rnadetector/ws/storage/app/tus_cache/" ] && mkdir -p "/rnadetector/ws/storage/app/tus_cache/"
[ ! -d "/rnadetector/ws/storage/app/logs/" ] && mkdir -p "/rnadetector/ws/storage/app/logs/"

create_data_dir
create_run_dir
initialize_mysql_database

chown -R www-data:staff "/rnadetector/ws"
chmod -R 777 "/rnadetector/ws/storage/"

[ "$DB_CREATED" = "true" ] && touch "${MYSQL_DATA_DIR}/ready"

[ -f /var/run/apache2/apache2.pid ] && rm -f /var/run/apache2/apache2.pid

echo "Starting supervisord"
exec supervisord -n
