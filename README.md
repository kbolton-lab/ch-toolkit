# CHIP Toolkit

A collection of utilities and tools used in the Clonal Hematopoiesis of Indeterminate Potential (CHIP) pipeline.  More details coming soon!

## Database Management

### Proposed CLI

    chip-variant-db init --db /path/to/sample.db --sample-name sample1
    chip-variant-db import-vcf --caller=lofreq --input-vcf=sample1.lofreq.vcf.gz /path/to/sample.db
    chip-variant-db dump-vcf --output-vcf=dump.vcf.gz /path/to/sample.db
    chip-variant-db merge-dbs --master-db=/path/to/master.db /path/to/sample.db
    chip-variant-db dump-stats --output-file /path/to/sample.db
    chip-variant-db update-db --input-vcf=sample1.pon_annotated.vcf.gz /path/to/sample.db

## Workflow Pipelines (Moonshot Project)

### Proposed CLI

#### Ingestion Pipeline

| Pipeline Attribute | Description |
| ------------------ | ----------- |
| **Goal:** | Normalize various sample inputs into either a standardized input format |
| **Main Input:**  | FASTQs or unaligned bams |
| **Main Output:** | An aligned and indexed bam. |

    chip-workflow ingestion \
        --input-fastqs fof-of-fastqs.csv \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml

#### Alignment and Variant Calling Pipeline

| Pipeline Attribute | Description |
| ------------------ | ----------- |
| **Goal:** | Align and variant call the indexed bam (per sample) |
| **Main Input:** | An indexed bam. |
| **Main Output:** | A database of variant calls per sample |

    chip-workflow align-vc \
        --input-bam indexed.bam \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml

#### Panel of Normal (PoN) Computations

| Pipeline Attribute | Description |
| ------------------ | ----------- |
| **Goal:** | Normalize sample inputs into either an indexed.bam _(or indexed.cram)_. |
| **Main Input:** | A database of sample variants. |
| **Main Output:** | A database of variant statistics |

_(This pipeline could work on a cohort of sample variants for improved performance.)_

    chip-workflow pon-computation \
        --input-variant-db /path/to/samples.db \
        --input-pon-db /path/to/pon.db \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml

#### Variant Annotations

| Pipeline Attribute | Description |
| ------------------ | ----------- |
| **Goal:** | Annotate a set of variants. |
| **Main Input:** | A database of variants from 1 or more samples. |
| **Main Output:** | An updated database of variants from 1 or more samples. |

_(This pipeline could work on a cohort of sample variants for improved performance.)_

    chip-workflow annotation \
        --input-db /path/to/samples.db \
        --input-pon-db /path/to/pon.db \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml

#### Machine Learning

| Pipeline Attribute | Description |
| ------------------ | ----------- |
| **Goal:** | Run candidate variants through a ML model |
| **Main Input:** | A database of variants from 1 or more samples. |
| **Main Output:** | An updated database of variants from 1 or more samples? (ML model dependent?) |

_(This pipeline could work on a cohort of sample variants for improved performance.)_

    chip-workflow ml-model \
        --input-db /path/to/samples.db \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml

#### Finalization

| Pipeline Attribute | Description |
| ------------------ | ----------- |
| **Goal:** | Assemble, Clean and Merge any disparate data or files or  Make Final Standarized Reports |
| **Main Input:** | A database of variants from 1 or more samples. |
| **Main Output:** | ??? |

_(This pipeline could work on a cohort of sample variants for improved performance.)_

    chip-workflow finalization \
        --input-db /path/to/samples.db \
        --output-directory /path/to/outdir \
        --configs /path/to/other-configs.yaml \
        --jobmgr /path/to/scheduler-configs.yaml
