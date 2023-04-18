import os, io, re
import pandas as pd

import chip.vdbtools.importers.vcf as vcf
import chip.vdbtools.importers.variants as variants
import chip.utils.logger as log
import chip.utils.database as db
from clint.textui import indent

def load_df_file_into_annotation(connection, df, table):
    sql = f"""
        INSERT INTO {table} SELECT df.*
        FROM df
        WHERE df.variant_id NOT IN (
            SELECT variant_id
            FROM {table} t
            WHERE t.variant_id IN (
                SELECT variant_id
                FROM df
            )
        )
    """
    log.logit(f"Starting to insert pandas dataframe into duckdb")
    connection.execute(sql)
    log.logit(f"Finished inserting pandas dataframe into duckdb")

def splitAndMax(row):
    res = {}
    # For each column, and value check if we need to split. If we do need to split, we need to find the maximum value
    # If not, then just return the float converted value (this eliminates the need for astype(float) earlier)
    for col, val in row.items():
        res[col] = max([float(re.sub(r'\.(?!\d)', '0', x)) for x in str(val).split(',')]) if ',' in str(val) else float(val)
    #res = {col : max([float(re.sub(r'\.(?!\d)', '0', x)) for x in str(val).split(',')]) if ',' in str(val) else float(val) for col, val in row.items()}
    return pd.Series(res)

#chr1:12828529:C:A variant_id: 71547 - For cases like this... we need to do some pre-processing because the original fixed_b38_exome.vcf.gz has duplicates
def fix_gnomADe(connection, debug):
    sql = f"""
            SELECT variant_id, gnomADe_AF, gnomADe_AF_AFR, gnomADe_AF_AMR, gnomADe_AF_ASJ, gnomADe_AF_EAS, gnomADe_AF_FIN, gnomADe_AF_NFE, gnomADe_AF_OTH, gnomADe_AF_SAS
            FROM vep
            WHERE gnomADe_AF LIKE '%,%' OR
                  gnomADe_AF_AFR LIKE '%,%' OR
                  gnomADe_AF_AMR LIKE '%,%' OR
                  gnomADe_AF_ASJ LIKE '%,%' OR
                  gnomADe_AF_EAS LIKE '%,%' OR
                  gnomADe_AF_FIN LIKE '%,%' OR
                  gnomADe_AF_NFE LIKE '%,%' OR
                  gnomADe_AF_OTH LIKE '%,%' OR
                  gnomADe_AF_SAS LIKE '%,%'
    """
    df = connection.execute(sql).df()
    total = len(df)
    log.logit(f"There were {total} variants that need to be fixed. Run using debug to see these variants.")
    if debug: df.to_csv('gnomADe_variants.csv', index=False)
    res = df.loc[:, df.columns != 'variant_id'].apply(lambda x: splitAndMax(x), axis=1)
    df = pd.concat([df['variant_id'], res], axis=1)
    sql = f"""
            UPDATE vep
            SET gnomADe_AF = df.gnomADe_AF,
                gnomADe_AF_AFR = df.gnomADe_AF_AFR,
                gnomADe_AF_AMR = df.gnomADe_AF_AMR,
                gnomADe_AF_ASJ = df.gnomADe_AF_ASJ,
                gnomADe_AF_EAS = df.gnomADe_AF_EAS,
                gnomADe_AF_FIN = df.gnomADe_AF_FIN,
                gnomADe_AF_NFE = df.gnomADe_AF_NFE,
                gnomADe_AF_OTH = df.gnomADe_AF_OTH,
                gnomADe_AF_SAS = df.gnomADe_AF_SAS
            FROM df
            WHERE vep.variant_id = df.variant_id
    """
    connection.execute(sql)

def annotateGnomad(df, debug):
    log.logit(f"Formatting gnomAD Information")
    df[df.filter(regex="^gnomAD[eg]*_.*").columns] = df[df.filter(regex="^gnomAD[eg]*_.*").columns].replace({"":0, ".":0, "-":0}).astype(float)
    df['max_gnomAD_AF_VEP'] = df.filter(regex=("^gnomAD_.*AF")).max(axis=1)
    df['max_gnomADe_AF_VEP'] = df.filter(regex=("^gnomADe_AF.")).max(axis=1)
    df['max_gnomADg_AF_VEP'] = df.filter(regex=("^gnomADg_AF.")).max(axis=1)
    df['max_pop_gnomAD_AF'] = df[['gnomAD_AF', 'gnomADe_AF', 'gnomADg_AF']].max(axis=1)
    return df

def prepareAnnotatePdData(df, debug):
    AminoAcids = {"Cys":"C", "Asp":"D", "Ser":"S", "Gln":"Q", "Lys":"K",
                  "Ile":"I", "Pro":"P", "Thr":"T", "Phe":"F", "Asn":"N",
                  "Gly":"G", "His":"H", "Leu":"L", "Arg":"R", "Trp":"W",
                  "Ala":"A", "Val":"V", "Glu":"E", "Tyr":"Y", "Met":"M",
                  "%3D":"=", "=":"="}
    return df
    #import pdb; pdb.set_trace()
    #df['HGVSp']
  #
  # final$AAchange <- gsub("(.*p\\.)(.*)", "\\2", final$HGVSp_VEP)
  # for (i in 1:length(AminoAcids)) {
  #   final$AAchange <- gsub(names(AminoAcids)[i], AminoAcids[i], final$AAchange)
  # }
  # final$gene_loci_p <- paste(final$SYMBOL_VEP,
  #   paste0(
  #     sapply(final$AAchange, function(x) str_split(x, "[0-9]+", n = 2)[[1]][1]),
  #     as.numeric(str_extract(final$AAchange, "\\d+"))
  #   ),
  #   sep = "_"
  # )
  # final$gene_loci_c <- paste(final$SYMBOL_VEP,
  #   gsub(".*:", "", final$HGVSc),
  #   sep = "_"
  # )
  # final$gene_loci_vep <- ifelse(is.na(final$gene_loci_p), final$gene_loci_c, final$gene_loci_p)
  # final$key <- with(final, paste(CHROM, POS, REF, ALT, sep = ":"))
  # final$gene_aachange <- with(final, paste(SYMBOL_VEP, AAchange, sep = "_"))
  # final$gene_cDNAchange <- paste(final$SYMBOL_VEP, gsub(".*:", "", final$HGVSc_VEP), sep = "_")
  #
  # ## ch_pd stuff
  # vars <- supportData[["vars"]]
  # vars$gene_aachange <- paste(vars$SYMBOL_VEP, vars$AAchange2, sep = "_")
  # vars$gene_cDNAchange <- paste(vars$SYMBOL_VEP, gsub(".*:", "", vars$HGVSc_VEP), sep = "_")
  # vars <- vars[vars$key %in% final$key | vars$gene_loci_vep %in% final$gene_loci_vep, ]
  #
  # dims <- dim(final)[[1]]
  # final <- sqldf("SELECT l.*, r.`n.loci.vep`, r.`source.totals.loci`
  #           FROM `final` as l
  #           LEFT JOIN `vars` as r
  #           on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep")
  # final <- final[!duplicated(final), ]
  #
  # ## make sure aachange exists as in doesn't end with an '_'; example: DNMT3A_ for splice
  # vars <- vars[!(grepl("_$", vars$gene_aachange) | grepl("_$", vars$gene_cDNAchange)), ]
  # final <- sqldf("SELECT l.*, r.`n.HGVSp`, r.`source.totals.p`
  #           FROM `final` as l
  #           LEFT JOIN `vars` as r
  #           on l.key = r.key OR (l.gene_aachange = r.gene_aachange)")
  # final <- final[!duplicated(final), ]
  # final <- sqldf("SELECT l.*, r.`n.HGVSc`, r.`source.totals.c`
  #           FROM `final` as l
  #           LEFT JOIN `vars` as r
  #           on l.key = r.key OR (l.gene_cDNAchange = r.gene_cDNAchange)")
  # final <- final[!duplicated(final), ]
  # paste0("dims match after sqldf: ", dim(final)[[1]] == dims)
  #
  # final$Mutect2_PON_2AT2_percent <- tidyr::replace_na(final$Mutect2_PON_2AT2_percent, 0)
  # final$Vardict_PON_2AT2_percent <- tidyr::replace_na(final$Vardict_PON_2AT2_percent, 0)

def preprocess(df, debug):
    log.logit(f"Performing some preprocessing to prepare for AnnotatePD")
    df = annotateGnomad(df, debug)
    df = prepareAnnotatePdData(df, debug)
    return df

def tsv_to_pd(vep, batch_number, debug):
    log.logit(f"Reading in the VEP VCF...")
    window = 5_000_000
    all_res = []
    header = None
    count = 0
    from itertools import islice
    with open(vep, 'rt') as f, indent(4, quote=' >'):
        while True:
            nextLines = list(islice(f, window))
            if not header:
                header = [l.strip('\n').split('\t') for l in nextLines if l.startswith('#Uploaded_variation')][0]
            nextLines = [l for l in nextLines if not l.startswith('#')]
            if not nextLines:
                break
            res = pd.read_csv(io.StringIO(''.join(nextLines)),
                                header=None,
                                sep='\t',
                                names=header,
                                low_memory=False).rename(columns={'#Uploaded_variation':'variant_id'})
            all_res.append(res)
            count += res.shape[0]
            cur = res.iloc[-1]['Location']
            log.logit(f"{count} variants loaded. Current: {cur}")
    all_res = pd.concat(all_res, ignore_index=True)
    total = len(all_res)
    log.logit(f"Finished reading in the VEP VCF: {total} variants.")
    all_res['batch'] = batch_number
    return total, all_res

def import_vep(variant_db, annotation_db, vep, batch_number, debug, clobber):
        log.logit(f"Adding VEP from batch: {batch_number} into {annotation_db}", color="green")
        annotation_connection = db.duckdb_connect_rw(annotation_db, clobber)
        extension = os.path.splitext(vep)[1]
        if extension == ".tsv":
            counts, df = tsv_to_pd(vep, batch_number, debug)
        else: # This should be deprecated (it's slower and worse)
            counts, df = vcf.vcf_to_pd(vep, "vep", batch_number, debug)
            variant_connection = db.duckdb_connect(variant_db)
            df = variants.insert_variant_id(df, variant_connection, debug)
            variant_connection.close()
        if annotation_connection.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='vep'").fetchone():
            log.logit(f"The VEP table already exists, so we can insert the information directly")
            load_df_file_into_annotation(annotation_connection, df, "vep")
        else:
            log.logit(f"This is the first time the VEP table is being referenced. Creating the Table.")
            annotation_connection.sql("CREATE TABLE IF NOT EXISTS vep AS SELECT * FROM df")
        annotation_connection.close()
        log.logit(f"Finished importing VEP information")
        log.logit(f"Variants Processed - Total: {counts}", color="green")
        log.logit(f"All Done!", color="green")

def annotate_pd(annotation_db, batch_number, debug, clobber):
    log.logit(f"Annotating variants from batch: {batch_number} in {annotation_db}", color="green")
    annotation_connection = db.duckdb_connect_rw(annotation_db, clobber)
    fix_gnomADe(annotation_connection, debug)
    log.logit(f"Grabbing Variants to perform AnnotatePD")
    df = annotation_connection.execute(f"SELECT * FROM vep WHERE batch = \'{batch_number}\'").df()
    df = preprocess(df, debug)
    #annotatePD
