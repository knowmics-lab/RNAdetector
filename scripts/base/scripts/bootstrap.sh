#!/usr/bin/env bash
set -e

[[ $DEBUG == true ]] && set -x

MYSQL_DATA_DIR="/rnadetector/ws/storage/app/database/"
RNADETECTOR_DB_ARCHIVE="/opt/database.tar.bz2"
MYSQL_USER="www-data"
MYSQL_GROUP="staff"
MYSQL_RUN_DIR="/var/run/mysqld"
MYSQL_LOG_FILE="/rnadetector/ws/storage/app/logs/mysqld.log"
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

create_log_dir() {
  if [ ! -f ${MYSQL_LOG_FILE} ]; then
    touch ${MYSQL_LOG_FILE}
  fi
  chmod -R 0775 ${MYSQL_LOG_FILE}
  chown -R ${MYSQL_USER}:${MYSQL_GROUP} ${MYSQL_LOG_FILE}
}

initialize_mysql_database() {
  # initialize MySQL data directory
  if [ ! -d "${MYSQL_DATA_DIR}/mysql" ]; then
    echo "Installing database..."
    if ! mysqld --initialize-insecure; then
      mysql_install_db >/dev/null 2>&1
    fi
  fi
  if [ ! -d "${MYSQL_DATA_DIR}/${DB_NAME}" ]; then
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

if [[ "$CLOUD_ENV" == "true" ]]; then
  echo "Starting RnaDetector Webservice in Cloud Mode"
  ## Disable mysql server startup
  [ ! -f /etc/supervisor/supervisord.conf.disabled ] && mv /etc/supervisor/supervisord.conf /etc/supervisor/supervisord.conf.disabled
  [ -f /etc/supervisor/supervisord-cloud.conf.disabled ] && mv /etc/supervisor/supervisord-cloud.conf.disabled /etc/supervisor/supervisord.conf
  [ -f /rnadetector/ws/.env ] && mv /rnadetector/ws/.env /rnadetector/ws/.env.disabled
  [ -f /rnadetector/ws/.env.cloud ] && mv /rnadetector/ws/.env.cloud /rnadetector/ws/.env
  [ -f /genkey.sh ] && rm /genkey.sh

  if [ ! -f "/rnadetector/ws/storage/app/.migrated" ] && [[ "$DEBUG" != "true" ]]; then
    if ! touch "/rnadetector/ws/storage/app/.migrated" &&
      php /rnadetector/ws/artisan migrate --seed --force &&
      php /rnadetector/ws/artisan first:boot; then
      rm "/rnadetector/ws/storage/app/.migrated"
    fi
  fi

else
  [ ! -f "/rnadetector/ws/storage/app/version_number" ] && [ -d "$MYSQL_DATA_DIR" ] && (rm -f $MYSQL_DATA_DIR/ib_logfile* || echo "Nothing to remove")
  create_data_dir
  create_run_dir
  create_log_dir
  initialize_mysql_database
fi

chown -R www-data:staff "/rnadetector/ws" &
chmod -R 777 "/rnadetector/ws/storage/" &

if [[ "$CLOUD_ENV" != "true" ]]; then
  [ "$DB_CREATED" = "true" ] && touch "${MYSQL_DATA_DIR}/ready"
fi

[ -f /var/run/apache2/apache2.pid ] && rm -f /var/run/apache2/apache2.pid

echo "Starting supervisord"
exec supervisord -n
