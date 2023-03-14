#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# INPUTS
BATCH=$1
DB=$2
CHROM=$3
REDIS_HOST=$4
REDIS_PORT=$5

# Command
VIRTUALENV=${SCRIPT_DIR}/../venv
BASE_CMD=${VIRTUALENV}/bin/chip-variant-db

# Redis

${BASE_CMD} ingest-variants \
    --db=${IN_VCF} \
    --redis-host=${REDIS_HOST} \
    --redis-port=${REDIS_PORT} \
    --batch-number=${BATCH} \
    --chromosome=${CHROM} \
    --clobber \
    --window-size 100000
