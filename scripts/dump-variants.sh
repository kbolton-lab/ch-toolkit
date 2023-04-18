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

${BASE_CMD} dump-variants \
    --vdb=${DB} \
    --header-type=${DUMMY_HEADER} \
    --batch-number=${BATCH} \
    --chromosome=${CHROM}
