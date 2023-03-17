#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# INPUTS
BATCH=$1
DB=$2
CHROM=$3
DUMMY_HEADER=$4

# Command
VIRTUALENV=${SCRIPT_DIR}/../venv_ic
BASE_CMD=${VIRTUALENV}/bin/chip-variant-db

# Redis

${BASE_CMD} dump-variants \
    --db=${DB} \
    --header=${DUMMY_HEADER}
    --batch-number=${BATCH}
    --chromosome=${CHROM}
