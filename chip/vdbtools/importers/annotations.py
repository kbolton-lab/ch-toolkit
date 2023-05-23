import os, io, re
import pandas as pd
import duckdb

import chip.vdbtools.importers.vcf as vcf
import chip.vdbtools.importers.variants as variants
import chip.utils.logger as log
import chip.utils.database as db
from clint.textui import indent

ANNOTATION_FILES="/storage1/fs1/bolton/Active/Protected/Annotation_Files/"
BOLTON_BICK_VARS="bick.bolton.vars3.txt"

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

#chr1:12828529:C:A variant_id: 71547 - For cases like this... we need to do some pre-processing because the original fixed_b38_exome.vcf.gz has duplicates
def fix_gnomADe(df, debug):
    for col in df.filter(regex='^gnomAD[eg]*_.*').columns:
        if debug: log.logit(f"Fixing: {col}")
        tmp = df[col].copy()
        tmp.replace(['-', '.'], None, inplace=True)
        tmp[~tmp.isnull()] = tmp[~tmp.isnull()].str.split(',', expand=True).replace('.', None).astype(float).max(axis=1)
        df[col] = tmp
    return df

def annotateGnomad(df, debug):
    log.logit(f"Formatting gnomAD Information...")
    #df[df.filter(regex="^gnomAD[eg]*_.*").columns] = df[df.filter(regex="^gnomAD[eg]*_.*").columns].replace({"":0, ".":0, "-":0}).astype(float)
    df['max_gnomAD_AF_VEP'] = df.filter(regex=("^gnomAD_.*AF")).max(axis=1)
    df['max_gnomADe_AF_VEP'] = df.filter(regex=("^gnomADe_AF.")).max(axis=1)
    df['max_gnomADg_AF_VEP'] = df.filter(regex=("^gnomADg_AF.")).max(axis=1)
    df['max_pop_gnomAD_AF'] = df[['gnomAD_AF', 'gnomADe_AF', 'gnomADg_AF']].max(axis=1)
    return df

def prepareAnnotatePdData(df, vars, debug):
    log.logit(f"Formatting VEP Information...")
    AminoAcids = {"Cys":"C", "Asp":"D", "Ser":"S", "Gln":"Q", "Lys":"K",
                  "Ile":"I", "Pro":"P", "Thr":"T", "Phe":"F", "Asn":"N",
                  "Gly":"G", "His":"H", "Leu":"L", "Arg":"R", "Trp":"W",
                  "Ala":"A", "Val":"V", "Glu":"E", "Tyr":"Y", "Met":"M",
                  "%3D":"=", "=":"="}

    df['AAchange'] = df['HGVSp'].str.extract(r'(.*p\.)(.*)')[1]
    df.replace({'AAchange': AminoAcids}, inplace=True, regex=True)
    df['loci_p'] = df["AAchange"].str.extract(r'(.[0-9]+)')
    df['gene_loci_p'] = df['SYMBOL']+"_"+df['loci_p']
    df['loci_c'] = df['HGVSc'].str.extract(r'(.*:)(.*)')[1]
    df['gene_loci_c'] = df['SYMBOL']+"_"+df['loci_c']
    df['gene_loci'] = df['gene_loci_p']
    df.loc[df['gene_loci_p'].isnull(),'gene_loci'] = df['gene_loci_c']
    df['gene_aachange'] = df['SYMBOL']+"_"+df['AAchange']
    df['gene_cDNAchange'] = df['SYMBOL']+"_"+df['loci_c']

    #BOLTON_BICK_VARS Formatting
    vars['gene_aachange'] = vars['SYMBOL_VEP']+"_"+vars['AAchange2']
    vars['gene_cDNAchange'] = vars['SYMBOL_VEP']+"_"+vars['HGVSc_VEP'].str.extract(r'(.*:)(.*)')[1]

    log.logit(f"Merging information from {BOLTON_BICK_VARS}.", color="yellow")
    with indent(4, quote=' >'):
        dims = len(df)
        df_tmp = df[(df['key'].isin(vars['key'])) | (df['gene_loci'].isin(vars['gene_loci_vep'].dropna())) | (df['gene_aachange'].isin(vars['gene_aachange'].dropna())) | (df['gene_cDNAchange'].isin(vars['gene_cDNAchange'].dropna()))]
        if len(df_tmp) > 0:
            df_tmp = df_tmp[['key', 'gene_loci', 'gene_aachange', 'gene_cDNAchange']]
            df_tmp['truncating'] = "not"
            df_tmp.loc[df['AAchange'].str.contains("Ter", na=False), 'truncating'] = "truncating"

            tmp = vars[(vars['key'].isin(df_tmp['key'])) | (vars['gene_loci_vep'].isin(df_tmp['gene_loci'].dropna()))]
            sql = f'''
                    SELECT l.*, r."n.loci.vep", r."source.totals.loci"
                    FROM df_tmp l
                    LEFT JOIN tmp r
                    ON l.key = r.key OR l.gene_loci = r.gene_loci_vep
            '''
            df_tmp = duckdb.sql(sql).df()
            df_tmp.drop_duplicates(inplace=True)
            variants = len(df_tmp)
            log.logit(f"Adding n.loci to {variants} variants.")

            sql = f'''
                    SELECT l.*, r."n.loci.truncating.vep", r."source.totals.loci.truncating"
                    FROM df_tmp l
                    LEFT JOIN tmp r
                    ON (l.key = r.key AND l.truncating = r.truncating) OR (l.gene_loci = r.gene_loci_vep AND l.truncating = r.truncating)
            '''
            df_tmp = duckdb.sql(sql).df()
            df_tmp.drop_duplicates(inplace=True)
            df_tmp.drop(['truncating'], axis=1, inplace=True)
            variants = len(df_tmp)
            log.logit(f"Adding n.loci to {variants} variants.")

            tmp = vars[(vars['key'].isin(df_tmp['key'])) | (vars['gene_aachange'].isin(df_tmp['gene_aachange'].dropna()))]
            sql = f'''
                    SELECT l.*, r."n.HGVSp", r."source.totals.p"
                    FROM df_tmp l
                    LEFT JOIN tmp r
                    ON l.key = r.key OR (l.gene_aachange = r.gene_aachange)
            '''
            df_tmp = duckdb.sql(sql).df()
            df_tmp.drop_duplicates(inplace=True)
            variants = len(df_tmp)
            log.logit(f"Adding n.loci to {variants} variants.")

            tmp = vars[(vars['key'].isin(df_tmp['key'])) | (vars['gene_cDNAchange'].isin(df_tmp['gene_cDNAchange'].dropna()))]
            sql = f'''
                    SELECT l.*, r."n.HGVSc", r."source.totals.c"
                    FROM df_tmp l LEFT JOIN tmp r
                    ON l.key = r.key OR (l.gene_cDNAchange = r.gene_cDNAchange)
            '''
            df_tmp = duckdb.sql(sql).df()
            df_tmp.drop_duplicates(inplace=True)
            variants = len(df_tmp)
            log.logit(f"Adding n.loci to {variants} variants.")
        else:
            #bolton_bick_vars=['n.loci.vep','source.totals.loci', 'n.loci.truncating.vep', 'source.totals.loci.truncating', 'n.HGVSp', 'source.totals.p', 'n.HGVSc', 'source.totals.c']
            #df_tmp.reindex(columns=[*df.columns.tolist(), *bolton_bick_vars])
            df_tmp['n.loci.vep'] = None
            df_tmp['source.totals.loci'] = None
            df_tmp['n.loci.truncating.vep'] = None
            df_tmp['source.totals.loci.truncating'] = None
            df_tmp['n.HGVSp'] = None
            df_tmp['source.totals.p'] = None
            df_tmp['n.HGVSc'] = None
            df_tmp['source.totals.c'] = None

        df_tmp = df_tmp[['key', 'n.loci.vep', 'source.totals.loci', 'n.loci.truncating.vep', 'source.totals.loci.truncating', 'n.HGVSp', 'source.totals.p', 'n.HGVSc', 'source.totals.c']]
        #df_tmp.drop(['gene_loci', 'gene_aachange', 'gene_cDNAchange'], axis=1, inplace=True)
        log.logit(f"Summarizing information from {BOLTON_BICK_VARS}.", color="yellow")
        df = pd.merge(df, df_tmp, on=['key'], how='left')

        if len(df) != dims: log.logit(f"ERROR: Something went wrong in the join. Dimensions don't match!", color="red")
        if len(df) == dims: log.logit(f"SUCCESS.", color="yellow")
    df.drop(['loci_p', 'loci_c', 'Location'], axis=1, inplace=True)
    return df

def preprocess(df, debug):
    log.logit(f"Performing some preprocessing to prepare for AnnotatePD")
    vars = pd.read_csv(ANNOTATION_FILES+BOLTON_BICK_VARS, sep='\t')
    df = annotateGnomad(df, debug)
    df = prepareAnnotatePdData(df, vars, debug)
    return df

def tsv_to_pd(tsv, batch_number, header, debug):
    log.logit(f"Reading in the VEP VCF...")
    res = pd.read_csv(io.StringIO(''.join(tsv)),
                        header=None,
                        sep='\t',
                        names=header,
                        low_memory=False).rename(columns={'#Uploaded_variation':'variant_id'})
    total = len(res)
    log.logit(f"Finished reading {total} variants in the VEP VCF.")
    res['batch'] = batch_number
    return total, res

def annotate_pd_to_pd(annotate_pd, batch_number, header, debug):
    log.logit(f"Reading in the AnnotatePD CSV...")
    window = 5_000_000
    all_res = []
    total = 0
    from itertools import islice
    with open(annotate_pd, 'rt') as f, indent(4, quote=' >'):
        header = f.readline().strip('\n').replace('\"', '').split(',')
        while True:
            nextLines = list(islice(f, window))
            nextLines = [l for l in nextLines]
            if not nextLines:
                break
            res = pd.read_csv(io.StringIO(''.join(nextLines)),
                                header=None,
                                names=header,
                                low_memory=False)
            res['batch'] = batch_number
            total += len(res)
            log.logit(f"{total} variants loaded.")
            all_res.append(res)
    all_res = pd.concat(all_res, ignore_index=True)
    total = len(all_res)
    log.logit(f"Finished reading in the AnnotatePD CSV: {total} variants.")
    res['batch'] = batch_number
    return total, res

def insert_vep(vep, annotation_connection, variant_connection, batch_number, debug):
    window = 5_000_000
    total = 0
    from itertools import islice
    with open(vep, 'rt') as f, indent(4, quote=' >'):
        header = [l.strip('\n').split('\t') for l in f if l.startswith('#Uploaded_variation')][0]
        f.seek(0) # Reset the reader to top
        while True:
            nextLines = list(islice(f, window))
            nextLines = [l for l in nextLines if not l.startswith('#')]
            if not nextLines:
                break
            count, df = tsv_to_pd(nextLines, batch_number, header, debug)
            total += count
            cur = df.iloc[-1]['Location']
            log.logit(f"{total} variants loaded.")
            df = fix_gnomADe(df, debug)
            df = variants.insert_variant_keys(df, variant_connection, debug)
            df = preprocess(df, debug)
            df.drop(['WildtypeProtein', 'FrameshiftSequence'], axis=1, inplace=True)
            if annotation_connection.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='vep'").fetchone():
                log.logit(f"The VEP table already exists, so we can insert the information directly")
                annotation_connection.execute("SET GLOBAL pandas_analyze_sample=0")
                load_df_file_into_annotation(annotation_connection, df, "vep")
            else:
                log.logit(f"This is the first time the VEP table is being referenced. Creating the Table.")
                # Sometimes if the PD has too many NULL it cannot figure out the type to cast so it fails. See: https://github.com/duckdb/duckdb/issues/6811
                annotation_connection.execute("SET GLOBAL pandas_analyze_sample=0")
                annotation_connection.sql("CREATE TABLE IF NOT EXISTS vep AS SELECT * FROM df")
    return total

def import_vep(annotation_db, variant_db, vep, batch_number, debug, clobber):
        log.logit(f"Adding VEP from batch: {batch_number} into {annotation_db}", color="green")
        annotation_connection = db.duckdb_connect_rw(annotation_db, clobber)
        variant_connection = db.duckdb_connect_ro(variant_db)
        counts = insert_vep(vep, annotation_connection, variant_connection, batch_number, debug)
        annotation_connection.close()
        log.logit(f"Finished importing VEP information")
        log.logit(f"Variants Processed - Total: {counts}", color="green")
        log.logit(f"All Done!", color="green")

def dump_variants_batch(annotation_db, batch_number, debug):
    log.logit(f"Dumping variants from batch: {batch_number} in {annotation_db} to a CSV file for AnnotatePD.", color="green")
    annotation_connection = db.duckdb_connect_ro(annotation_db)
    log.logit(f"Grabbing Variants to perform AnnotatePD")
    sql = f'''
            SELECT variant_id, key, Consequence, SYMBOL, EXON, AAchange, HGVSc, HGVSp, \"n.HGVSc\", \"n.HGVSp\"
            FROM vep
            WHERE batch = {batch_number} AND SYMBOL != \'-\' AND Consequence NOT LIKE 'intron_variant%'
    '''
    df = annotation_connection.execute(sql).df()
    df[['CHROM', 'POS', 'REF', 'ALT']] = df['key'].str.split(':', expand=True)
    df.rename({'SYMBOL':'SYMBOL_VEP',
               'HGVSp':'HGVSp_VEP',
               'HGVSc':'HGVSc_VEP',
               'Consequence':'Consequence_VEP',
               'EXON':'EXON_VEP'}, axis=1, inplace=True)
    total = len(df)
    log.logit(f"{total} variants grabbed from {annotation_db} for {batch_number}")
    annotation_connection.close()
    df.to_csv(f"batch-{batch_number}-forAnnotatePD.csv", index=False)
    log.logit(f"Finished dumping variants into CSV file")
    log.logit(f"All Done!", color="green")

def process_annotate_pd(df, debug):
    df.drop(['Consequence_VEP', 'SYMBOL_VEP', 'EXON_VEP', 'AAchange', 'HGVSc_VEP', 'HGVSp_VEP', 'n.HGVSc', 'n.HGVSp', 'CHROM', 'POS', 'REF', 'ALT'], axis=1, inplace=True)
    return df

def import_annotate_pd(annotation_db, annotate_pd, batch_number, debug):
    log.logit(f"Adding AnnotatePD Information from batch: {batch_number} into {annotation_db}", color="green")
    annotation_connection = db.duckdb_connect_rw(annotation_db, False)
    counts, df = annotate_pd_to_pd(annotate_pd, annotation_connection, batch_number, debug)
    df = process_annotate_pd(df, debug)
    if annotation_connection.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='pd'").fetchone():
        log.logit(f"The AnnotatePD table already exists, so we can insert the information directly")
        annotation_connection.execute("SET GLOBAL pandas_analyze_sample=0")
        load_df_file_into_annotation(annotation_connection, df, "pd")
    else:
        log.logit(f"This is the first time the AnnotatePD table is being referenced. Creating the Table.")
        # Sometimes if the PD has too many NULL it cannot figure out the type to cast so it fails. See: https://github.com/duckdb/duckdb/issues/6811
        annotation_connection.execute("SET GLOBAL pandas_analyze_sample=0")
        annotation_connection.sql("CREATE TABLE IF NOT EXISTS pd AS SELECT * FROM df")
    log.logit(f"Finished importing AnnotatePD information")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")
