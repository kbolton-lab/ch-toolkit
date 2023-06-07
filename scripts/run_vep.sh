#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

function log {
    local timestamp=$(date +"%Y-%m-%d %T")
    echo "---> [ ${timestamp} ] $@" >&2
}

# LSF Job Parameters
DOCKER_IMAGE='docker(kboltonlab/ic_vep:latest)'
MEMORY=32GB
RUSAGE="rusage[mem=${MEMORY}]"
SELECT="select[mem>${MEMORY}]"
SPAN='span[hosts=1]'
COMPUTE_GROUP='compute-bolton'
QUEUE='general'
JOB_GROUP='/chani/multi'

export LSF_DOCKER_VOLUMES="/home/$USER:/home/$USER /storage1/fs1/bolton/Active:/storage1/fs1/bolton/Active /scratch1/fs1/bolton:/scratch1/fs1/bolton /storage1/fs1/bga/Active:/storage1/fs1/bga/Active"

# INPUTS
#COMBINED_VCF=$1
#DB=$2
#CHROM=$3
#DUMMY_HEADER=$4
ROOT_LOG_DIR=/scratch1/fs1/bolton/chani/chip-toolkit/run-vep
CHROMS=($(seq 1 22))
CHROMS+=( "X" "Y" )

echo "===> Running VEP <==="
i=1
for chrom in ${CHROMS[@]}; do
    job_name="VEP.chr${chrom}.vcf.gz)"
    logfile=${ROOT_LOG_DIR}/VEP.chr${chrom}.%J.log
    vcf=/storage1/fs1/bolton/Active/Projects/chip-toolkit/data/derived/3.1-dump-variants/chr${chrom}.vcf.gz
    outfile=chr${chrom}_VEP_annotated.tsv
    log "[ ${i} ] Processing chromosome: ${chrom}"
    log "VCF File: ${vcf}"
    log "Log Directory: ${ROOT_LOG_DIR}"
    log "LSF Job Name: ${job_name}"
    bsub \
        -J ${job_name} \
        -M ${MEMORY} \
        -R "${RUSAGE} ${SELECT} ${SPAN}" \
        -G ${COMPUTE_GROUP} \
        -g ${JOB_GROUP} \
        -q ${QUEUE} \
        -a "${DOCKER_IMAGE}" \
        -o ${logfile} \
        /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
            --fork 4 -i ${vcf} --tab -o ${outfile} \
            --offline --cache --buffer_size 1000 \
            --symbol --transcript_version --assembly GRCh38 --cache_version 104 --species homo_sapiens --merged --use_given_ref \
            --synonyms /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/chromAlias.ensembl.txt \
            --pick --pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length \
            --dir /storage1/fs1/bolton/Active/Projects/mocha/UKBB/ukbb_calls/pvcf/vep_zip \
            --fasta ${HG38_REF_DH} \
            --af_gnomad \
            --plugin Frameshift --plugin Wildtype \
            --plugin SpliceAI,snv=/storage1/fs1/bolton/Active/Protected/Data/hg38/vcf/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/storage1/fs1/bolton/Active/Protected/Data/hg38/vcf/spliceai_scores.raw.indel.hg38.vcf.gz \
            --everything \
            --check_existing --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/fixed_b38_exome.vcf.gz,gnomADe,vcf,exact,1,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
            --check_existing --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/custom_clinvar_vcf/v20181028/custom.vcf.gz,clinvar,vcf,exact,1,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES \
            --check_existing --custom /storage1/fs1/bolton/Active/Protected/Data/hg38/vcf/gnomad.genomes.r3.0.sites.full_trimmed_info.af_only.vcf.gz,gnomADg,vcf,exact,1,AF,AF_ami,AF_oth,AF_afr,AF_sas,AF_asj,AF_fin,AF_amr,AF_nfe,AF_eas \
            --force_overwrite --no_stats
    ((i += 1))
done
