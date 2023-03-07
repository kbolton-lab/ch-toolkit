#!/bin/bash

DOCKER_IMAGE='docker(indraniel/bolton-db-toolkit:v1)'
BATCH=1
MEMORY='30GB'


LSF_DOCKER_VOLUMES="/home/$USER:/home/$USER /storage1/fs1/bolton/Active:/storage1/fs1/bolton/Active /scratch1/fs1/bolton:/scratch1/fs1/bolton" \
  bsub \
  -Is \
  -M ${MEMORY} \
  -R "'rusage[mem=${MEMORY},tmp=30] select[mem>32GB && tmp>30] span[hosts=1]'" \
  -G compute-bolton \
  -q general-interactive \
  -a "'${DOCKER_IMAGE}'" \
  /bin/bash -l
