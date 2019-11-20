#!/usr/bin/env bash

service nginx start
service php7.2-fpm start
service supervisor start

exec "$@"
