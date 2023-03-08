#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# INPUTS
BATCH=$1
IN_VCF=$2
REDIS_HOST=$3
REDIS_PORT=$4

# Command
VIRTUALENV=${SCRIPT_DIR}/../venv
BASE_CMD=${VIRTUALENV}/bin/chip-variant-db

# Redis

${BASE_CMD} register-variants \
    --input-vcf=${IN_VCF} \
    --redis-host=${REDIS_HOST} \
    --redis-port=${REDIS_PORT} \
    --batch-number=${BATCH}
