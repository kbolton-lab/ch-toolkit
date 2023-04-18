#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

function log {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "---> [ ${timestamp} ] $@" >&2
}

BATCH=1
DUMMY_HEADER="dummy"

# LSF Job Parameters
DOCKER_IMAGE='docker(indraniel/bolton-db-toolkit:v1)'
MEMORY=2GB
RUSAGE="rusage[mem=${MEMORY}]"
SELECT="select[mem>2GB]"
SPAN='span[hosts=1]'
COMPUTE_GROUP='compute-bolton'
QUEUE='general'
JOB_GROUP='/chani/multi'
SCRIPT=${SCRIPT_DIR}/dump-variants.sh

export LSF_DOCKER_VOLUMES="/home/$USER:/home/$USER /storage1/fs1/bolton/Active:/storage1/fs1/bolton/Active /scratch1/fs1/bolton:/scratch1/fs1/bolton"

# input gathering
CHROMS=($(seq 1 22))
CHROMS+=( "X" "Y" )


ROOT_LOG_DIR=/scratch1/fs1/bolton/chani/chip-toolkit/dump-variants

echo "===> Submitting chromosomes <==="

i=1
for chrom in ${CHROMS[@]}; do
    job_name="chr${chrom}.dump.variants.batch.${BATCH}"
    logfile=${ROOT_LOG_DIR}/batch${BATCH}.chr${chrom}.%J.log
    db=/storage1/fs1/bolton/Active/Projects/chip-toolkit/variants.db
    log "[ ${i} ] Processing chromosome: ${chrom}"
    log "Database File: ${db}"
    log "Log Directory: ${ROOT_LOG_DIR}"
    log "LSF Job Name: ${job_name}"
    set -o xtrace;
    if (( i > 1 )); then
        bsub \
        -w "ended(${jobid})" \
        -J ${job_name} \
        -M ${MEMORY} \
        -R "${RUSAGE} ${SELECT} ${SPAN}" \
        -G ${COMPUTE_GROUP} \
        -g ${JOB_GROUP} \
        -q ${QUEUE} \
        -a "${DOCKER_IMAGE}" \
        -o ${logfile} \
        /bin/bash ${SCRIPT} ${BATCH} ${db} "chr${chrom}" ${DUMMY_HEADER}
    else
        bsub \
        -J ${job_name} \
        -M ${MEMORY} \
        -R "${RUSAGE} ${SELECT} ${SPAN}" \
        -G ${COMPUTE_GROUP} \
        -g ${JOB_GROUP} \
        -q ${QUEUE} \
        -a "${DOCKER_IMAGE}" \
        -o ${logfile} \
        /bin/bash ${SCRIPT} ${BATCH} ${db} "chr${chrom}" ${DUMMY_HEADER}
    fi
    jobid=$(bjobs -o "jobid" -J ${job_name} | awk 'NR==2{print $1}')
    set +o xtrace;
    ((i += 1))
done
