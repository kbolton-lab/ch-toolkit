# CHIP Toolkit

A collection of utilities and tools used in the CHIP pipeline.  More details coming soon!

## Database Management

### Proposed CLI

    chip-db init-db /path/to/sample.db # initialize a variant database
    chip-db import-vcf --caller=lofreq --input-vcf=sample1.lofreq.vcf.gz /path/to/sample.db
    chip-db dump-vcf --output-vcf=dump.vcf.gz /path/to/sample.db
    chip-db merge-dbs --master-db=/path/to/master.db /path/to/sample.db
    chip-db dump-stats --output-file /path/to/sample.db
