import os, glob
import vcfpy
import duckdb
import pandas as pd

import chip.vdbtools.importers.vcf as vcf
import chip.vdbtools.importers.variants as variants
import chip.utils.logger as log
import chip.utils.database as db
import chip.utils.fisher_exact_test as fisher_test
from clint.textui import indent

def ensure_mutect_tbl(connection):
    log.logit("Ensuring or creating the mutect table")
    sql = '''
        DROP TYPE IF EXISTS mutect_filter_type;
        CREATE TYPE mutect_filter_type AS ENUM (
            'PASS',
            'FAIL',
            'base_qual',
            'clustered_events',
            'contamination',
            'duplicate',
            'fragment',
            'germline',
            'haplotype',
            'low_allele_frac',
            'map_qual',
            'multiallelic',
            'n_ratio',
            'normal_artifact',
            'orientation',
            'panel_of_normals',
            'position',
            'possible_numt',
            'slippage',
            'strand_bias',
            'strict_strand',
            'weak_evidence'
        );
        CREATE TABLE IF NOT EXISTS mutect(
            sample_id                           integer NOT NULL,
            variant_id                          integer NOT NULL,
            version                             varchar(50),
            mutect_filter                       mutect_filter_type[],
            info_as_filterstatus                varchar[],       /* e.g. AS_FilterStatus=strand_bias|strand_bias */
            info_as_sb_table                    varchar(20),     /* e.g. AS_SB_TABLE=5379,6627|1,27|1,174 */
            info_dp                             integer,
            info_ecnt                           integer,
            info_mbq_ref                        integer,
            info_mbq_alt                        integer,
            info_mfrl_ref                       integer,
            info_mfrl_alt                       integer,
            info_mmq_ref                        integer,
            info_mmq_alt                        integer,
            info_mpos                           integer,
            info_popaf                          decimal(10,1),
            info_roq                            decimal(10,1),
            info_rpa_ref                        integer,
            info_rpa_alt                        integer,
            info_ru                             varchar(10),
            info_str                            boolean,        /* 0 Flag */
            info_strq                           integer,
            info_tlod                           decimal(10,5),
            format_af                           decimal(10,5),
            format_dp                           integer,
            format_ref_count                    integer,
            format_alt_count                    integer,
            format_ref_f1r2                     integer,
            format_alt_f1r2                     integer,
            format_ref_f2r1                     integer,
            format_alt_f2r1                     integer,
            format_gt                           varchar(3),
            format_ref_fwd                      integer,
            format_ref_rev                      integer,
            format_alt_fwd                      integer,
            format_alt_rev                      integer,
            pon_2at2_percent                    boolean,
            pon_nat2_percent                    integer,
            pon_max_vaf                         decimal(10,5),
            fisher_p_value                      decimal(22,20),
            PRIMARY KEY( sample_id, variant_id )
        )
    '''
    connection.execute(sql)

def ensure_vardict_tbl(connection):
    log.logit("Ensuring or creating the vardict table")

    sql = '''
        DROP TYPE IF EXISTS vardict_filter_type;
        CREATE TYPE vardict_filter_type AS ENUM (
            'PASS',
            'q22.5',
            'Q10',
            'p8',
            'SN1.5',
            'Bias',
            'pSTD',
            'd3',
            'v2',
            'min_af',
            'MSI12',
            'NM5.25',
            'InGap',
            'InIns',
            'Cluster0bp',
            'LongMSI',
            'AMPBIAS',
            'BCBIO'
        );
        CREATE TABLE IF NOT EXISTS vardict(
            sample_id               integer NOT NULL,
            variant_id              integer NOT NULL,
            version                 varchar(10),
            vardict_filter          vardict_filter_type[],
            info_type               varchar(5),
            info_dp                 integer,
            info_vd                 integer,
            info_af                 decimal(10,5),
            info_bias               varchar(20),      /* 11671:3913 */
            info_refbias            varchar(20),
            info_varbias            varchar(20),
            info_pmean              decimal(10,5),
            info_pstd               decimal(10,5),
            info_qual               decimal(10,5),
            info_qstd               decimal(10,5),
            info_sbf                decimal(10,5),
            info_oddratio           decimal(10,5),
            info_mq                 decimal(10,5),
            info_sn                 decimal(10,5),
            info_hiaf               decimal(10,5),
            info_adjaf              decimal(10,5),
            info_shift3             integer,
            info_msi                decimal(10,5),
            info_msilen             decimal(10,5),
            info_nm                 decimal(10,5),
            info_lseq               varchar(20),
            info_rseq               varchar(20),
            info_hicnt              integer,
            info_hicov              integer,
            info_splitread          integer,
            info_spanpair           integer,
            info_duprate            decimal(10,5),
            format_ref_count        integer,
            format_alt_count        integer,
            format_gt               varchar(3),
            format_dp               integer,
            format_vd               integer,
            format_af               decimal(10,5),
            format_ref_fwd          integer,
            format_ref_rev          integer,
            format_alt_fwd          integer,
            format_alt_rev          integer,
            pon_2at2_percent        boolean,
            pon_nat2_percent        integer,
            pon_max_vaf             decimal(10,5),
            fisher_p_value          decimal(22,20),
            PRIMARY KEY( sample_id, variant_id )
        )
    '''
    connection.execute(sql)

def drop_indexes(connection):
    log.logit("Dropping existing indexes on the caller table")
    sql = "DROP INDEX IF EXISTS sample_variant_idx"
    connection.execute(sql)

def create_indexes(connection, caller):
    if caller == "mutect":
        log.logit("Creating new indexes on the mutect table")
        log.logit("Generating the sample_variant_idx")
        sql = """
            CREATE UNIQUE INDEX sample_variant_idx
            ON mutect ( sample_id, variant_id )
        """
    elif caller == "vardict":
        log.logit("Creating new indexes on the vardict table")
        log.logit("Generating the sample_variant_idx")
        sql = """
            CREATE UNIQUE INDEX sample_variant_idx
            ON vardict ( sample_id, variant_id )
        """
        connection.execute(sql)

def setup_caller_tbl(connection, caller):
    if caller == "mutect":
        log.logit("Preparing the mutect database file")
        ensure_mutect_tbl(connection)
    elif caller == "vardict":
        log.logit("Preparing the vardict database file")
        ensure_vardict_tbl(connection)
    drop_indexes(connection)
    return connection

def process_mutect(df, variant_connection, debug):
    log.logit(f"Formatting Mutect Dataframe...")
    df['version'] = '2.2'
    df = df.drop(columns=['CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL'])
    df = df.rename({'FILTER': 'mutect_filter',
                    'info_pon_2at2_percent':'pon_2at2_percent',
                    'info_pon_nat2_percent':'pon_nat2_percent',
                    'info_pon_max_vaf':'pon_max_vaf'}, axis=1)
    df.loc[df['info_as_filterstatus'].isnull(), 'info_as_filterstatus'] = 'Multiallelic'
    df["info_as_filterstatus"] = df["info_as_filterstatus"].apply(lambda x: [x])
    df.loc[df['info_str'].isnull(), 'info_str'] = False
    df[['info_mbq_ref', 'info_mbq_alt']] = df['info_mbq'].str.split(',', expand=True)
    df[['info_mfrl_ref', 'info_mfrl_alt']] = df['info_mfrl'].str.split(',', expand=True)
    df[['info_mmq_ref', 'info_mmq_alt']] = df['info_mmq'].str.split(',', expand=True)
    df[['info_rpa_ref', 'info_rpa_alt']] = df['info_rpa'].str.split(',', expand=True)
    df[['format_ref_count', 'format_alt_count']] = df['format_ad'].str.split(',', expand=True)
    df[['format_ref_f1r2', 'format_alt_f1r2']] = df['format_f1r2'].str.split(',', expand=True)
    df[['format_ref_f2r1', 'format_alt_f2r1']] = df['format_f2r1'].str.split(',', expand=True)
    df[['format_ref_fwd', 'format_ref_rev', 'format_alt_fwd', 'format_alt_rev']] = df['format_sb'].str.split(',', expand=True)
    df['fisher_p_value'] = None
    log.logit(f"Adding Variant IDs to the Mutect Variants")
    df = variants.insert_variant_id(df, variant_connection, debug)
    df = df[['sample_id', 'variant_id', 'version', 'mutect_filter', 'info_as_filterstatus', 'info_as_sb_table', 'info_dp', 'info_ecnt', 'info_mbq_ref', 'info_mbq_alt', 'info_mfrl_ref', 'info_mfrl_alt', 'info_mmq_ref', 'info_mmq_alt', 'info_mpos', 'info_popaf', 'info_roq', 'info_rpa_ref', 'info_rpa_alt', 'info_ru', 'info_str', 'info_strq', 'info_tlod', 'format_af', 'format_dp', 'format_ref_count', 'format_alt_count', 'format_ref_f1r2', 'format_alt_f1r2', 'format_ref_f2r1', 'format_alt_f2r1', 'format_gt', 'format_ref_fwd', 'format_ref_rev', 'format_alt_fwd', 'format_alt_rev', 'pon_2at2_percent', 'pon_nat2_percent', 'pon_max_vaf', 'fisher_p_value']]
    return df

def insert_mutect_caller(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug):
    log.logit(f"Registering mutect variants from: {input_vcf} in batch: {batch_number}", color="green")
    sample_name = os.path.basename(input_vcf).split('.')[1]
    log.logit(f"Grabbing sample_id from {sample_db} for sample: {sample_name}")
    sample_id = db.duckdb_connect(sample_db).execute(f"SELECT sample_id FROM samples WHERE sample_name = '{sample_name}';").fetchone()[0]
    counts, df = vcf.vcf_to_pd(input_vcf, "caller", batch_number, debug)
    df['sample_id'] = sample_id
    variant_connection = db.duckdb_connect(variant_db)
    df = process_mutect(df, variant_connection, debug)
    variant_connection.close()
    caller_connection = db.duckdb_connect_rw(db_path, clobber)
    setup_caller_tbl(caller_connection, "mutect")
    vcf.duckdb_load_df_file(caller_connection, df, "mutect")
    #create_indexes(caller_connection, "mutect")
    caller_connection.close()
    annotate_fisher_test(variant_db, db_path, "mutect", batch_number, debug)
    log.logit(f"Finished inserting mutect variants")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def process_vardict(df, variant_connection, debug):
    log.logit(f"Formatting Vardict Dataframe...")
    df['version'] = 'v1.8.2'
    df = df.drop(columns=['CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL'])
    df = df.rename({'FILTER': 'vardict_filter',
                    'info_pon_2at2_percent':'pon_2at2_percent',
                    'info_pon_nat2_percent':'pon_nat2_percent',
                    'info_pon_max_vaf':'pon_max_vaf'}, axis=1)
    df[['format_ref_count', 'format_alt_count']] = df['format_ad'].str.split(',', expand=True)
    df[['format_ref_fwd', 'format_ref_rev']] = df['format_rd'].str.split(',', expand=True)
    df[['format_alt_fwd', 'format_alt_rev']] = df['format_ald'].str.split(',', expand=True)
    df['fisher_p_value'] = None
    log.logit(f"Adding Variant IDs to the Vardict Variants")
    df = variants.insert_variant_id(df, variant_connection, debug)
    df = df[['sample_id', 'variant_id', 'version', 'vardict_filter', 'info_type', 'info_dp', 'info_vd', 'info_af', 'info_bias', 'info_refbias', 'info_varbias', 'info_pmean', 'info_pstd', 'info_qual', 'info_qstd', 'info_sbf', 'info_oddratio', 'info_mq', 'info_sn', 'info_hiaf', 'info_adjaf', 'info_shift3', 'info_msi', 'info_msilen', 'info_nm', 'info_lseq', 'info_rseq', 'info_hicnt', 'info_hicov', 'info_splitread', 'info_spanpair', 'info_duprate', 'format_ref_count', 'format_alt_count', 'format_gt', 'format_dp', 'format_vd', 'format_af', 'format_ref_fwd', 'format_ref_rev', 'format_alt_fwd', 'format_alt_rev', 'pon_2at2_percent', 'pon_nat2_percent', 'pon_max_vaf', 'fisher_p_value']]
    return df

def insert_vardict_caller(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug):
    log.logit(f"Registering vardict variants from: {input_vcf} in batch: {batch_number}", color="green")
    sample_name = os.path.basename(input_vcf).split('.')[1]
    log.logit(f"Grabbing sample_id from {sample_db} for sample: {sample_name}")
    sample_id = db.duckdb_connect(sample_db).execute(f"SELECT sample_id FROM samples WHERE sample_name = '{sample_name}';").fetchone()[0]
    counts, df = vcf.vcf_to_pd(input_vcf, "caller", batch_number, debug)
    df['sample_id'] = sample_id
    variant_connection = db.duckdb_connect(variant_db)
    df = process_vardict(df, variant_connection, debug)
    variant_connection.close()
    caller_connection = db.duckdb_connect_rw(db_path, clobber)
    setup_caller_tbl(caller_connection, "vardict")
    vcf.duckdb_load_df_file(caller_connection, df, "vardict")
    #create_indexes(caller_connection, "vardict")
    caller_connection.close()
    annotate_fisher_test(variant_db, db_path, "vardict", batch_number, debug)
    log.logit(f"Finished inserting vardict variants")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def merge_caller_tables(db_path, connection, batch_number, caller, debug):
    log.logit(f'Merging Sample Callers from {db_path}')
    with indent(4, quote=' >'):
        for i, file in enumerate(glob.glob(db_path + "/" + "*.db")):
            sample_name = os.path.basename(file).split('.')[1]
            log.logit(f"Merging: {file}")
            connection.execute(f"ATTACH \'{file}\' as sample_{i}")
            connection.sql(f"INSERT INTO {caller} SELECT * FROM sample_{i}.{caller}")
    log.logit(f"Finished merging all tables from: {db_path}")

def ingest_caller_batch(db_path, caller_db, caller, batch_number, debug, clobber):
    log.logit(f"Ingesting variants from batch: {batch_number} into {caller_db}", color="green")
    connection = db.duckdb_connect_rw(caller_db, clobber)
    setup_caller_tbl(connection, caller)
    merge_caller_tables(db_path, connection, batch_number, caller, debug)
    connection.close()
    log.logit(f"Finished ingesting variants")
    log.logit(f"All Done!", color="green")

def annotate_fisher_test(variant_db, caller_db, caller, batch_number, debug):
    log.logit(f"Performing the Fisher's Exact Test on the variants inside {caller_db} for batch: {batch_number}")
    caller = "mutect" if caller.lower() == "mutect" else "vardict"
    caller_connection = db.duckdb_connect_rw(caller_db, False)
    log.logit(f"Finding all variants within {caller_db} that does not have the fisher's exact test p-value calcualted...")
    caller_connection.execute(f"ATTACH '{variant_db}' as temp_variants")
    sql = f'''
        SELECT c.variant_id, c.sample_id, v.PoN_RefDepth, v.PoN_AltDepth, c.format_ref_fwd, c.format_ref_rev, c.format_alt_fwd, c.format_alt_rev,
        FROM {caller} c LEFT JOIN temp_variants.variants v
        ON c.variant_id = v.variant_id
        WHERE c.fisher_p_value is NULL AND PoN_RefDepth is NOT NULL AND PoN_AltDepth is NOT NULL
    '''
    df = caller_connection.execute(sql).df()
    length = len(df)
    log.logit(f"There were {length} variants without fisher test p-value within {caller_db}")
    log.logit(f"Calculating Fisher Exact Test for all variants inside {caller_db}")
    df = fisher_test.pvalue_df(df)
    log.logit(f"Updating {caller_db} with the fisher's exact test p-values")
    sql = f"""
        UPDATE {caller} as c
        SET fisher_p_value = df.pvalue
        FROM df
        WHERE c.variant_id = df.variant_id AND c.sample_id = df.sample_id
    """
    caller_connection.sql(sql)
    caller_connection.close()
    log.logit(f"Finished updating fisher test p-values inside {caller_db}")
    log.logit(f"Done!", color = "green")
