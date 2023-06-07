#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# INPUTS
IN_VCF=$1
BATCH=$2
DB=$3

# Command
VIRTUALENV=${SCRIPT_DIR}/../venv_ic
BASE_CMD=${VIRTUALENV}/bin/chip-variant-db

${BASE_CMD} register-sample-variants \
    --input-vcf=${IN_VCF} \
    --batch-number=${BATCH} \
    --db=${DB} \
    --clobber
