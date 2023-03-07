#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

function log {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "---> [ ${timestamp} ] $@" >&2
}

# LSF Job Parameters
DOCKER_IMAGE='docker(indraniel/bolton-db-toolkit:v1)'
MEMORY=30GB
TMP_SPACE=30
RUSAGE="rusage[mem=${MEMORY},tmp=${TMP_SPACE}]"
SELECT="select[mem>32GB && tmp>${TMP_SPACE}]"
SPAN='span[hosts=1]'
COMPUTE_GROUP='compute-bolton'
QUEUE='general-interactive'
SCRIPT=${SCRIPT_DIR}/register-variants.sh

export LSF_DOCKER_VOLUMES="/home/$USER:/home/$USER /storage1/fs1/bolton/Active:/storage1/fs1/bolton/Active /scratch1/fs1/bolton:/scratch1/fs1/bolton"

# Registration Command Parameters
BATCH=1
ATTEMPT=1

REDIS_HOST=compute1-exec-130.ris.wustl.edu
REIDS_PORT=8999

VCF_ROOT_DIR=/storage1/fs1/bolton/Active/Projects/chip-toolkit/data/TERRA_Data
MUTECT_VCFS=$(find ${VCF_ROOT_DIR} -name "mutect.*.vcf.gz" -print)
NUM_MUTECT=${#MUTECT_VCFS[@]}
VARDICT_VCFS=$(find ${VCF_ROOT_DIR} -name "vardict.*.vcf.gz" -print)
NUM_VARDICT=${#VARDICT_VCFS[@]}

ROOT_LOG_DIR=/scratch1/fs1/bolton/idas/variant-registration-logs

"===> Submitting mutect vcfs <==="

i=0
for vcf in ${MUTECT_VCFS[@]}; do
    sample=$(basename $(dirname ${vcf}))
    job_name="${sample}.mutect.batch.${BATCH}.${ATTEMPT}"
    log_dir=${ROOT_LOG_DIR}/mutect/${sample}
    mkdir -p ${log_dir}
    logfile="${log_dir}/%J.log"
    log "[ ${i} | ${NUM_MUTECT}] Processing VCF: ${vcf}"
    log "Sample: ${sample}"
    log "Log Directory: ${log_dir}"
    log "LSF Job Name: ${job_name}"
    LSF_DOCKER_VOLUMES=
      bsub \
      -Is \
      -J ${job_name} \
      -M ${MEMORY} \
      -R "${RUSAGE} ${SELECT} ${SPAN}" \
      -G ${COMPUTE_GROUP} \
      -q ${QUEUE} \
      -a "${DOCKER_IMAGE}" \
      -o ${logfile} \
      ${SCRIPT} ${BATCH} ${VCF} ${REDIS_HOST} ${REDIS_PORT}
     ((i += 1))
done

"===> Submitting vardict vcfs <==="

i=0
for vcf in ${VARDICT_VCFS[@]}; do
    sample=$(basename $(dirname ${vcf}))
    job_name="${sample}.vardict.batch.${BATCH}.${ATTEMPT}"
    log_dir=${ROOT_LOG_DIR}/mutect/${sample}
    mkdir -p ${log_dir}
    logfile="${log_dir}/%J.log"
    log "[ ${i} | ${NUM_VARDICT}] Processing VCF: ${vcf}"
    log "Sample: ${sample}"
    log "Log Directory: ${log_dir}"
    log "LSF Job Name: ${job_name}"
    LSF_DOCKER_VOLUMES=
      bsub \
      -Is \
      -J ${job_name} \
      -M ${MEMORY} \
      -R "${RUSAGE} ${SELECT} ${SPAN}" \
      -G ${COMPUTE_GROUP} \
      -q ${QUEUE} \
      -a "${DOCKER_IMAGE}" \
      -o ${logfile} \
      ${SCRIPT} ${BATCH} ${VCF} ${REDIS_HOST} ${REDIS_PORT}
      ((i += 1))
done