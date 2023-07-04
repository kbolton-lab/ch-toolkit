# CH Toolkit

The CH Toolkit is a collection of utilities and tools created by Irenaeus Chan (chani@wustl.edu) and Indraniel for the purpose of handling and maintaing variants called by the [ArCH Pipeline](https://github.com/kbolton-lab/ArCH/tree/main).

## Sample Alignment &amp; Variant Calling Pipeline

Variant Calling is performed using two separate Variant Callers: GATK's Mutect2 and AstraZeneca's VarDict

To improve on cost and runtime performance, only genomic regions where Mutect2 successfully called and passed variants are used as input BED windows to VarDict.

Details behind this WDL workflow can be found here under the [ArCH WGS Variant Calling Pipeline](https://github.com/kbolton-lab/ArCH/blob/main/WDL/WGS/Subworkflows/variant_calling.wdl).

The subsequent variants called by both callers are then sanitized, normalized, and filtered for common Germline variants before being used as inputs for this toolkit.

## Variant Annotation

To consolidate as much information as possible prior to filtering for CH variants. Several annotation steps must be performed including:
* VEP - The effect of variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions.
* Panel of Normal (PoN) Pileup - Estimation of the statistical noise within a given region or location from sequencing
* AnnotatePD (RScript) - The classification of CH variants through prior large-scale sequencing projects and population data including OncoKB, COSMIC, gnomAD, etc.

## Variant Processing Pipeline

While this tool is meant to be run as a stand-alone infrastructure, for the purposes of automation, a [WDL workflow](https://github.com/kbolton-lab/ArCH/blob/main/WDL/WGS/ArCCH_WGS.wdl) has been written with [inputs](https://github.com/kbolton-lab/ArCH/blob/main/WDL/WGS/ArCCH_WGS.json) being simply the VCF files produced from Mutect2 and VarDict, sample phenotype information, and the local database filepath for the storage and maintainence of the variant database.

### CLI

```
Usage: ch-toolkit [OPTIONS] COMMAND [ARGS]...

  A collection of db related tools for handling sample data.

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  calculate-fishers-test  Updates the variants inside Mutect or Vardict tables
                          with p-value from Fisher's Exact Test
  dump-annotations        dumps all variant annotations inside duckdb into a
                          CSV file
  dump-ch                 Outputs CH Variants from Database
  dump-variants           dumps all variants inside duckdb into a VCF file
  import-annotate-pd      annotates variants with their pathogenicity
  import-pon-pileup       updates variants inside duckdb with PoN pileup
                          information
  import-sample-variants  Register the variants for a VCF file into a variant
                          database
  import-sample-vcf       import a vcf file into sample variant database
  import-samples          Loads a CSV containing samples into samples database
  import-vep              updates variants inside duckdb with VEP information
  merge-batch-variants    Combines all sample variant databases into a single
                          database
  merge-batch-vcf         Combines all sample vcfs databases into a single
                          database
```

### Workflow Diagram

![CH WGS Annotation Pipeline](docs/images/ch-pipeline-wgs.png)

## Functions

### Importing Samples
| import-samples ||
|-----------|-----------------------------------------------------------------------------------------|
|**Goal:**  | Create the ***samples.db*** database which will contain the information for the samples |
|**Input:** | A CSV file containing the information relevant to the samples being processed           |
|**Output:**| ***samples.db*** database                                                               |

```
  ch-toolkit import-samples \
    --samples washu-cad-1.csv \
    --sdb database/samples.db \
    --batch-number 1
```

### Importing Variants for a Single Sample
| import-sample-variants ||
|-----------|-----------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Create individual ***sample.variant.db*** database for each sample that will contain the variant specific information |
|**Input:** | Sample VCF file containing the information about the variants                                                         |
|**Output:**| ***sample.variant.db*** database                                                                                      |

```
  ch-toolkit import-sample-variants \
    --input-vcf mutect.sample_name.vcf.gz \
    --vdb variant_databases/mutect.sample_name.db \
    --batch-number 1

  ch-toolkit import-sample-variants \
    --input-vcf vardict.sample_name.vcf.gz \
    --vdb variant_databases/vardict.sample_name.db \
    --batch-number 1
```

### Create a Single Centralized Variants Database from all Individual Sample Variants Databases
| merge-batch-variants ||
|-----------|---------------------------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Create a single ***variants.db*** that will contain ALL unique variants found within the individual ***sample.variant.db*** databases |
|**Input:** | Path where the ***sample.varaint.db*** database files are located                                                                     |
|**Output:**| ***variants.db*** database                                                                                                            |

```
  ch-toolkit merge-batch-variants \
    --db-path variant_databases/ \
    --vdb database/variants.db \
    --batch-number 1
```

### Output Variants to VCF
| dump-variants ||
|-----------|--------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Convert variants within ***variants.db*** into a VCF file for use in various annotation like VEP, PoN Pileup, etc. |
|**Input:** | ***variants.db*** database containing all unique variants that need to be converted into a VCF file                |
|**Output:**| A VCF file containing all unique variants from ***variants.db***                                                   |

```
  ch-toolkit dump-variants \
    --vdb database/variants.db \
    --header-type simple \
    --batch-number 1
```

### Importing VCF/Caller Specific Information for a Single Sample
| import-sample-vcf ||
|-----------|----------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Create individual ***sample.caller.db*** databases for each sample that will contain the caller specific information |
|**Input:** | Sample VCF file containing the information from the specific caller                                                  |
|**Output:**| ***sample.caller.db*** database                                                                                      |

```
  ch-toolkit import-sample-vcf \
    --caller mutect \
    --input-vcf mutect.sample_name.vcf.gz \
    --cdb mutect_databases/mutect.sample_name.db \
    --batch-number 1

  ch-toolkit import-sample-vcf \
    --caller vardict \
    --input-vcf vardict.sample_name.vcf.gz \
    --cdb vardict_databases/vardict.sample_name.db \
    --batch-number 1
```

### Create a Single Centralized Mutect and Vardict Caller Database from all Individual Sample Caller Databases
| merge-batch-vcf ||
|-----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Create a single ***mutect.db*** and ***vardict.db*** database that will contain ALL caller specific information found within the individual ***sample.caller.db*** databases |
|**Input:** | Path where the ***sample.caller.db*** database files are located                                                                                                             |
|**Output:**| ***mutect.db*** database or ***vardict.db*** database                                                                                                                        |

```
  ch-toolkit merge-batch-vcf \
    --db-path mutect_databases/ \
    --cdb database/mutect.db \
    --vdb database/variants.db \
    --sdb database/samples.db \
    --caller mutect \
    --batch-number 1	

  ch-toolkit merge-batch-vcf \
    --db-path mutect_databases/ \
    --cdb database/vardict.db \
    --vdb database/variants.db \
    --sdb database/samples.db \
    --caller vardict \
    --batch-number 1
```

### After Annotating Variants Using VEP, Import VEP Information
| import-vep ||
|-----------|--------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Import annotation information produced from annotating variants from the **dump-variants** step using VEP          |
|**Input:** | The resulting TSV file produced from running VEP using (`--tab`) mode                                              |
|**Output:**| ***annotations.db*** database                                                                                      |

```
  ch-toolkit import-vep \
    --adb database/annotations.db \
    --vdb database/variants.db \
    --vep VEP_annotated.tsv \
    --batch-number 1
```

### Leveraging the Annotations from VEP, output Variants that are Potentially Putative Drivers
| dump-annotations ||
|-----------|--------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Export variants that are potentially putative drivers to be annotated by the custom AnnotatePD RScript             |
|**Input:** | ***annotations.db*** database containing information about the variants that need to be converted into a CSV file  |
|**Output:**| A CSV file containing all unique variants from ***annotations.db*** that will be annotated using AnnotatePD        |

```
  ch-toolkit dump-annotations \
    --adb database/annotations.db \
    --batch-number 1
```

### After Annotating Variants Using AnnotatePD, Import AnnotatePD Information
| import-annotate-pd ||
|-----------|---------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Import annotation information produced from annotating variants from the **dump-annotations** step using AnnotatePD |
|**Input:** | The resulting CSV file produced from running AnnotatePD                                                             |
|**Output:**| ***annotations.db*** database                                                                                       |

```
  ch-toolkit import-annotate-pd \
    --adb database/annotations.db \
    --pd annotatePD_results.csv \
    --batch-number 1
```

### After Performing the [PoN Pileup Workflow](https://github.com/kbolton-lab/ArCH/blob/main/WDL/WGS/Subworkflows/PoN.wdl), Import Pileup Information
| import-pon-pileup ||
|-----------|--------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Import pileup information produced from running the [PoN Pileup Workflow](https://github.com/kbolton-lab/ArCH/blob/main/WDL/WGS/Subworkflows/PoN.wdl) on the variants from the **dump-variants** step          |
|**Input:** | The resulting PoN Pileup VCF file produced from the PoN Pileup Workflow                                            |
|**Output:**| ***pileup.db*** database                                                                                           |

```
  ch-toolkit import-pon-pileup \
    --vdb database/variants.db \
    --pdb database/pileup.db \
    --pon-pileup pon_pileup.vcf.gz \
    --batch-number 1
```

### Calculate the Signal-to-Noise Ratio for All Variants
| calculate-fishers-test ||
|-----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Perform a Fisher's Exact Test comparing the proportion of variant alleles to reference alleles to the proportion of the same variant alleles to reference alleles found within the PoN Samples |
|**Input:** | The pileup information within the ***pileups.db*** database along with variant information found within ***mutect.db*** or ***vardict.db*** databases                                          |
|**Output:**| ***mutect.db*** and ***vardict.db*** databases updated with a p-value indiciating the signifance of the detected variant signal relative to the expected noise for the same given location     |

```
  ch-toolkit calculate-fishers-test \
    --pdb database/pileup.db \
    --cdb database/mutect.db \
    --caller mutect \
    --batch-number 1

  ch-toolkit calculate-fishers-test \
    --pdb database/pileup.db \
    --cdb database/vardict.db \
    --caller vardict \
    --batch-number 1
```

### Filter and Identify Putative Driver Variants
| dump-ch ||
|-----------|-------------------------------------------------------------------------------------------------------------------------------|
|**Goal:**  | Process through all information currently stored in the databases and detect CH variants with pathogenic support              |
|**Input:** | ***mutect.db***, ***vardict.db***, and ***annotations.db*** databases                                                         |
|**Output:**| A CSV file containing all of the variants and all relevant information pertaining to said variants predicted to be pathogenic |

```
  ch-toolkit dump-ch \
    --mcdb database/mutect.db \
    --vcdb database/vardict.db \
    --adb database/annotations.db
```
