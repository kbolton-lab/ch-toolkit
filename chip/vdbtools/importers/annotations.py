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
TRUNCATING="BB.truncating.more.than.1.tsv"
SEGEMENTAL_DUPLICATIONS="dup.grch38.bed.gz"
SIMPLE_REPEATS="simpleRepeat.bed"
REPEAT_MASKER="repeatMaskerJoinedCurrent.bed"
MUT2_BICK="topmed.n2.mutation.c.p.txt"
MUT2_KELLY="kelly.n2.mutation.c.p.txt"
MATCHES2="matches.2.c.p.txt"
GENE_LIST="oncoKB_CGC_pd_table_disparity_KB_BW.csv"
ONCO_KB_AVAILABLE_GENE="/home/tran.n/DucTran/UkbbVar/Data/oncoKbCancerGeneList.tsv"

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

# DEPRECATED
def splitAndMax(row):
    res = {}
    # For each column, and value check if we need to split. If we do need to split, we need to find the maximum value
    # If not, then just return the float converted value (this eliminates the need for astype(float) earlier)
    for col, val in row.items():
        res[col] = max([float(re.sub(r'\.(?!\d)', '0', x)) for x in str(val).split(',')]) if ',' in str(val) else float(val)
    #res = {col : max([float(re.sub(r'\.(?!\d)', '0', x)) for x in str(val).split(',')]) if ',' in str(val) else float(val) for col, val in row.items()}
    return pd.Series(res)

#chr1:12828529:C:A variant_id: 71547 - For cases like this... we need to do some pre-processing because the original fixed_b38_exome.vcf.gz has duplicates
# DEPRECATED
def fix_gnomADe_(connection, debug):
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

    # Sometimes if the PD has too many NULL it cannot figure out the type to cast so it fails. See: https://github.com/duckdb/duckdb/issues/6811
    #temp_connection.execute("SET GLOBAL pandas_analyze_sample=0")
    log.logit(f"Merging information from {BOLTON_BICK_VARS}.", color="yellow")
    with indent(4, quote=' >'):
        dims = len(df)
        tmp = df[(df['key'].isin(vars['key'])) | (df['gene_loci'].isin(vars['gene_loci_vep'].dropna()))]
        variants = len(tmp)
        if variants > 0:
            sql = f'''
                    SELECT l.key, l.gene_loci, r."n.loci.vep", r."source.totals.loci"
                    FROM tmp l
                    LEFT JOIN vars r
                    ON l.key = r.key OR l.gene_loci = r.gene_loci_vep
            '''
            log.logit(f"Adding n.loci to {variants} variants.")
            tmp = duckdb.sql(sql).df()
            #tmp = temp_connection.execute(sql).df()
            tmp = tmp[['key', 'gene_loci', 'n.loci.vep', 'source.totals.loci']]
            tmp.drop_duplicates(inplace=True)
            df = pd.merge(df, tmp, on=['key', 'gene_loci'], how='left')
        else:
            df['n.loci.vep'] = None
            df['source.totals.loci'] = None

        df['truncating'] = "not"
        df.loc[df['AAchange'].str.contains("Ter", na=False), 'truncating'] = "truncating"
        tmp = df[(df['key'].isin(vars['key'])) | (df['gene_loci'].isin(vars['gene_loci_vep'].dropna()))]
        variants = len(tmp)
        if variants > 0:
            sql = f'''
                    SELECT l.key, l.gene_loci, l.truncating, r."n.loci.truncating.vep", r."source.totals.loci.truncating"
                    FROM tmp l
                    LEFT JOIN vars r
                    ON (l.key = r.key AND l.truncating = r.truncating) OR (l.gene_loci = r.gene_loci_vep AND l.truncating = r.truncating)
            '''
            log.logit(f"Adding n.loci.truncating.vep to {variants} variants.")
            #tmp = temp_connection.execute(sql).df()
            tmp = duckdb.sql(sql).df()
            tmp = tmp[['key', 'gene_loci', 'truncating', 'n.loci.truncating.vep', 'source.totals.loci.truncating']]
            tmp.drop_duplicates(inplace=True)
            df = pd.merge(df, tmp, on=['key', 'gene_loci', 'truncating'], how='left')
        else:
            df['n.loci.truncating.vep'] = None
            df['source.totals.loci.truncating'] = None
        df.drop(['truncating'], axis=1, inplace=True)

        tmp = df[(df['key'].isin(vars['key'])) | (df['gene_aachange'].isin(vars['gene_aachange'].dropna()))]
        variants = len(tmp)
        if variants > 0:
            sql = f'''
                    SELECT l.key, l.gene_aachange, r."n.HGVSp", r."source.totals.p"
                    FROM tmp l
                    LEFT JOIN vars r
                    ON l.key = r.key OR (l.gene_aachange = r.gene_aachange)
            '''
            log.logit(f"Adding n.HGVSp to {variants} variants.")
            #tmp = temp_connection.execute(sql).df()
            tmp = duckdb.sql(sql).df()
            tmp = tmp[['key', 'gene_aachange', 'n.HGVSp', 'source.totals.p']]
            tmp.drop_duplicates(inplace=True)
            df = pd.merge(df, tmp, on=['key', 'gene_aachange'], how='left')
        else:
            df['n.HGVSp'] = None
            df['source.totals.p'] = None

        tmp = df[(df['key'].isin(vars['key'])) | (df['gene_cDNAchange'].isin(vars['gene_cDNAchange'].dropna()))]
        variants = len(tmp)
        if variants > 0:
            sql = f'''
                    SELECT l.key, l.gene_cDNAchange, r."n.HGVSc", r."source.totals.c"
                    FROM tmp l
                    LEFT JOIN vars r
                    ON l.key = r.key OR (l.gene_cDNAchange = r.gene_cDNAchange)
            '''
            log.logit(f"Adding n.HGVSc to {variants} variants.")
            #tmp = temp_connection.execute(sql).df()
            tmp = duckdb.sql(sql).df()
            tmp = tmp[['key', 'gene_cDNAchange', 'n.HGVSc', 'source.totals.c']]
            tmp.drop_duplicates(inplace=True)
            df = pd.merge(df, tmp, on=['key', 'gene_cDNAchange'], how='left')
        else:
            df['n.HGVSc'] = None
            df['source.totals.c'] = None

        #temp_connection.close()
        #os.unlink(tmp_path)
        if len(df) != dims: log.logit(f"ERROR: Something went wrong in the join. Dimensions don't match!", color="red")
        if len(df) == dims: log.logit(f"SUCCESS.", color="green")
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
            log.logit(f"{total} variants loaded. Current: {cur}")
            df = fix_gnomADe(df, debug)
            df = variants.insert_variant_keys(df, variant_connection, debug)
            df = preprocess(df, debug)
            if annotation_connection.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='vep'").fetchone():
                log.logit(f"The VEP table already exists, so we can insert the information directly")
                load_df_file_into_annotation(annotation_connection, df, "vep")
            else:
                log.logit(f"This is the first time the VEP table is being referenced. Creating the Table.")
                annotation_connection.execute("SET GLOBAL pandas_analyze_sample=0")
                annotation_connection.sql("CREATE TABLE IF NOT EXISTS vep AS SELECT * FROM df")
            #fix_gnomADe_(annotation_connection, debug)
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
    log.logit(f"Annotating variants from batch: {batch_number} in {annotation_db} with putative driver information", color="green")
    annotation_connection = db.duckdb_connect_ro(annotation_db)
    log.logit(f"Grabbing Variants to perform AnnotatePD")
    sql = f'''
            SELECT variant_id, key, Consequence, SYMBOL, EXON, AAchange, HGVSc, HGVSp, \"n.HGVSc\", \"n.HGVSp\"
            FROM vep
            WHERE batch = {batch_number} AND SYMBOL != \'-\' AND Consequence NOT LIKE 'intron_variant%'
    '''
    df = annotation_connection.execute(sql).df()
    df[['CHROM', 'POS', 'REF', 'ALT']] = df['key'].str.split(':', expand=True)
    df.rename({'variant_id':'SAMPLE',
               'SYMBOL':'SYMBOL_VEP',
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
