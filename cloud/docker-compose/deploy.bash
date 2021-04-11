#!/usr/bin/env bash

WHITE='\033[1;37m'
NC='\033[0m'

# Define environment variables...
export APP_SERVICE=${APP_SERVICE:-"rnadetector"}

# Ensure we  are within the docker-compose environment
if [ ! -f ./docker-compose.yml ]; then
  echo -e "${WHITE}You should run this utility from the same directory as your docker-compose.yml file.${NC}" >&2
  exit 1
fi

# Ensure that Docker is running...
if ! docker info >/dev/null 2>&1; then
  echo -e "${WHITE}Docker is not running.${NC}" >&2
  exit 1
fi

if [ $# -gt 0 ]; then
  # Source deploy.conf file
  if [ -f ./deploy.conf ]; then
    source ./deploy.conf
    export APP_URL
    export APP_KEY
    export APP_PORT
    export DB_PASSWORD
  fi
  if [ "$1" == "mysql" ]; then
    shift 1
    if ! docker-compose up -d mysql; then
      echo -e "${WHITE}Unable to start mysql container${NC}" >&2
      exit 1
    fi
  elif [ "$1" == "app" ]; then
    shift 1
    if ! docker-compose up -d "$APP_SERVICE"; then
      echo -e "${WHITE}Unable to start app container${NC}" >&2
      exit 1
    fi
  elif [ "$1" == "exec" ]; then
    shift 1
    docker-compose exec "$APP_SERVICE" "$@"
  else
    docker-compose "$@"
  fi
else
  docker-compose ps
fi
