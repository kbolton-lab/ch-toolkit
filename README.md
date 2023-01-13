# CHIP Toolkit

A collection of utilities and tools used in the CHIP pipeline.  More details coming soon!

## Database Management

### Proposed CLI

    chip-variant-db init --db /path/to/sample.db --sample-name sample1
    chip-variant-db import-vcf --caller=lofreq --input-vcf=sample1.lofreq.vcf.gz /path/to/sample.db
    chip-variant-db dump-vcf --output-vcf=dump.vcf.gz /path/to/sample.db
    chip-variant-db merge-dbs --master-db=/path/to/master.db /path/to/sample.db
    chip-variant-db dump-stats --output-file /path/to/sample.db

## Workflow Pipelines (Moonshot Project)

### Proposed CLI

    chip-workflow ingestion \
        --input-fastqs fof-of-fastqs.csv \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml

    chip-workflow align-vc \
        --input-bams fof-of-bams.csv \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml

    chip-workflow pon-computation \
        --input-db /path/to/samples.db \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml

    chip-workflow annotation \
        --input-db /path/to/samples.db \
        --input-pon-db /path/to/pon.db \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml

    chip-workflow ml-model \
        --input-db /path/to/samples.db \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml

    chip-workflow finalization \
        --input-db /path/to/samples.db \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml
