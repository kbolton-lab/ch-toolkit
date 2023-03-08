#!/bin/bash

# interactive mode
LSF_DOCKER_PORTS='8999:8999' \
LSF_DOCKER_VOLUMES="/home/$USER:/home/$USER /storage1/fs1/bolton/Active:/storage1/fs1/bolton/Active /scratch1/fs1/bolton:/scratch1/fs1/bolton" \
  bsub \
  -Is \
  -M 30GB \
  -R 'rusage[mem=30GB,tmp=30] select[mem>32GB && tmp>30] span[hosts=1]' \
  -R 'select[port8999=1]' \
  -G compute-bolton \
  -q general-interactive \
  -a 'docker(indraniel/bolton-db-toolkit:v1)' \
  /bin/bash -l

# # batch mode
# LSF_DOCKER_PORTS='8999:8999' \
# LSF_DOCKER_VOLUMES="/home/$USER:/home/$USER /storage1/fs1/bolton/Active:/storage1/fs1/bolton/Active /scratch1/fs1/bolton:/scratch1/fs1/bolton" \
#   bsub \
#   -M 30GB \
#   -R 'rusage[mem=30GB,tmp=30] select[mem>32GB && tmp>30] span[hosts=1]' \
#   -R 'select[port8999=1]' \
#   -G compute-bolton \
#   -q general-interactive \
#   -o /scratch1/fs1/bolton/redis/batch-server-logs
#   -a 'docker(indraniel/bolton-db-toolkit:v1)' \
#   /opt/bolton-lab/redis-7.0.8/bin/redis-server /scratch1/fs1/bolton/redis/redis.conf --port 8999 --dbfilename 2023-03-07.rdb
