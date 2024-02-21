import os, glob, shutil, math
import duckdb
import pandas as pd
import multiprocessing as mp

import ch.vdbtools.handlers.vcf as vcf
import ch.vdbtools.handlers.variants as variants
import ch.vdbtools.handlers.samples as samples
import ch.utils.logger as log
import ch.utils.database as db
import ch.utils.fisher_exact_test as fisher_test
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
            sample_name                         varchar(50) NOT NULL,
            key                                 VARCHAR NOT NULL,
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
            sample_id                           integer,
            variant_id                          BIGINT,
            batch                               integer
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
            sample_name             varchar(50) NOT NULL,
            key                     VARCHAR NOT NULL,
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
            info_oddratio           decimal(30,5),
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
            sample_id               integer,
            variant_id              BIGINT,
            batch                   integer,
        )
    '''
    connection.execute(sql)

def setup_caller_tbl(connection, caller):
    if caller == "mutect":
        log.logit("Preparing the mutect database file")
        ensure_mutect_tbl(connection)
    elif caller == "vardict":
        log.logit("Preparing the vardict database file")
        ensure_vardict_tbl(connection)
    return connection

def process_mutect(df, debug):
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
    df['sample_id'] = None
    df['variant_id'] = None
    df = df[['sample_name', 'key', 'version', 'mutect_filter', 'info_as_filterstatus', 'info_as_sb_table', 'info_dp', 'info_ecnt', 'info_mbq_ref', 'info_mbq_alt', 'info_mfrl_ref', 'info_mfrl_alt', 'info_mmq_ref', 'info_mmq_alt', 'info_mpos', 'info_popaf', 'info_roq', 'info_rpa_ref', 'info_rpa_alt', 'info_ru', 'info_str', 'info_strq', 'info_tlod', 'format_af', 'format_dp', 'format_ref_count', 'format_alt_count', 'format_ref_f1r2', 'format_alt_f1r2', 'format_ref_f2r1', 'format_alt_f2r1', 'format_gt', 'format_ref_fwd', 'format_ref_rev', 'format_alt_fwd', 'format_alt_rev', 'pon_2at2_percent', 'pon_nat2_percent', 'pon_max_vaf', 'fisher_p_value', 'sample_id', 'variant_id', 'batch']]
    log.logit(f"Removing any duplicate variants that may have occured during the merging process")
    df = df.drop_duplicates(subset=['key', 'sample_name'], keep='first')
    return df

def insert_mutect_caller(db_path, input_vcf, batch_number, clobber, debug):
    log.logit(f"Registering mutect variants from: {input_vcf} in batch: {batch_number}", color="green")
    counts, df = vcf.vcf_to_pd(input_vcf, "caller", batch_number, debug)
    if 'info_sample' in df.columns:
        df['sample_name'] = df['info_sample']
    else:
        if input_vcf.endswith('.vcf'):
            df['sample_name'] = '.'.join(os.path.basename(input_vcf).split('.')[1:-1])
        elif input_vcf.endswith('.vcf.gz'):
            df['sample_name'] = '.'.join(os.path.basename(input_vcf).split('.')[1:-2])
    df = process_mutect(df, debug)
    caller_connection = db.duckdb_connect_rw(db_path, clobber)
    setup_caller_tbl(caller_connection, "mutect")
    vcf.duckdb_load_df_file(caller_connection, df, "mutect")
    #! Write out parquet file
    #parquetPath = db_path.replace(".db", ".parquet")
    #cmd = f"COPY mutect TO '{parquetPath}' (FORMAT 'parquet')"
    #caller_connection.execute(cmd)
    caller_connection.close()
    log.logit(f"Finished inserting mutect variants")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def process_vardict(df, debug):
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
    df['sample_id'] = None
    df['variant_id'] = None
    df = df[['sample_name', 'key', 'version', 'vardict_filter', 'info_type', 'info_dp', 'info_vd', 'info_af', 'info_bias', 'info_refbias', 'info_varbias', 'info_pmean', 'info_pstd', 'info_qual', 'info_qstd', 'info_sbf', 'info_oddratio', 'info_mq', 'info_sn', 'info_hiaf', 'info_adjaf', 'info_shift3', 'info_msi', 'info_msilen', 'info_nm', 'info_lseq', 'info_rseq', 'info_hicnt', 'info_hicov', 'info_splitread', 'info_spanpair', 'info_duprate', 'format_ref_count', 'format_alt_count', 'format_gt', 'format_dp', 'format_vd', 'format_af', 'format_ref_fwd', 'format_ref_rev', 'format_alt_fwd', 'format_alt_rev', 'pon_2at2_percent', 'pon_nat2_percent', 'pon_max_vaf', 'fisher_p_value', 'sample_id', 'variant_id', 'batch']]
    log.logit(f"Removing any duplicate variants that may have occured during the merging process")
    df = df.drop_duplicates(subset=['key', 'sample_name'], keep='first')
    return df

def insert_vardict_caller(db_path, input_vcf, batch_number, clobber, debug):
    log.logit(f"Registering vardict variants from: {input_vcf} in batch: {batch_number}", color="green")
    counts, df = vcf.vcf_to_pd(input_vcf, "caller", batch_number, debug)
    if 'info_sample' in df.columns:
        df['sample_name'] = df['info_sample']
    else:
        if input_vcf.endswith('.vcf'):
            df['sample_name'] = '.'.join(os.path.basename(input_vcf).split('.')[1:-1])
        elif input_vcf.endswith('.vcf.gz'):
            df['sample_name'] = '.'.join(os.path.basename(input_vcf).split('.')[1:-2])
    df = process_vardict(df, debug)
    caller_connection = db.duckdb_connect_rw(db_path, clobber)
    setup_caller_tbl(caller_connection, "vardict")
    vcf.duckdb_load_df_file(caller_connection, df, "vardict")
    #! Write out parquet file
    #parquetPath = db_path.replace(".db", ".parquet")
    #cmd = f"COPY vardict TO '{parquetPath}' (FORMAT 'parquet')"
    #caller_connection.execute(cmd)
    caller_connection.close()
    log.logit(f"Finished inserting vardict variants")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

#! DEPRECATED
# This is the merging using *.db - For merging using *.parquet, see below
def merge_caller_tables_(db_path, caller_connection, variant_db, sample_db, batch_number, caller, debug):
    log.logit(f'Merging Sample Callers from {db_path}')

    # Getting all sample names from the caller table
    sql = f"""
        SELECT DISTINCT sample_id
        FROM {caller}
    """
    caller_connection.execute(sql)
    sample_ids = caller_connection.df()["sample_id"].tolist()

    with indent(4, quote=' >'):
        for i, file in enumerate(glob.glob(db_path + "/" + "*.db")):
            sample_name = os.path.basename(file).split('.')[1]            
            with indent(4, quote=' >'):
                sample_caller_connection = db.duckdb_connect_rw(file, False)
                log.logit(f"Adding Sample IDs to {caller} database")
                samples.insert_sample_id_into_db(sample_caller_connection, caller, sample_db, debug)
                log.logit(f"Adding Variant IDs to the {caller} variants")
                variants.insert_variant_id_into_db(sample_caller_connection, caller, variant_db, debug)
                sample_caller_connection.close()
            log.logit(f"Merging: {i} - {file}")
            caller_connection.execute(f"ATTACH \'{file}\' as sample_{i} (READ_ONLY)")
            if debug: log.logit(f"Attached {file} as sample_{i}")

            # Getting Sample Name from Sample_VCF
            sql = f"""
                SELECT sample_id
                FROM sample_{i}.{caller}
                LIMIT 1
            """
            caller_connection.execute(sql)
            sample_id = caller_connection.df().iloc[0, 0]

            # If the sample exists, then we need to check for repeating variant_id and sample_id
            # If it doesn't exist, we can just import the information directly without checking.
            if sample_id in sample_ids:
                log.logit(f"Sample {sample_id} already exists in {caller}, using safe merge")
                sql = f"""
                    INSERT INTO {caller} SELECT s.*
                    FROM sample_{i}.{caller} s
                    WHERE s.variant_id NOT IN (
                        SELECT variant_id
                        FROM {caller} c
                        WHERE c.sample_id = '{sample_id}'
                    )
                """
            else:
                log.logit(f"Sample {sample_id} does not exist in {caller}, using unsafe merge")
                sql = f"""
                    INSERT INTO {caller} SELECT *
                    FROM sample_{i}.{caller}
                """
                sample_ids.append(sample_id)
            if debug: log.logit(f"Executing: {sql}")
            caller_connection.sql(sql)
            caller_connection.execute(f"DETACH sample_{i}")
            if debug: log.logit(f"SQL Complete")
    log.logit(f"Finished merging all tables from: {db_path}")

# Parallelize adding the SampleID and VariantID and copying out the Parquet File
def prepare_parquet_files(index, file, caller, sample_ids, sample_db, variant_db, check_folder, no_check_folder):
    check = False
    log.logit(f"Processing: {index} - {file}")
    sample_caller_connection = db.duckdb_connect_rw(file, False)
    #sample_caller_connection.execute(f"ALTER TABLE {caller} ALTER COLUMN variant_id TYPE BIGINT")
    log.logit(f"Adding Sample IDs to {caller} database")
    samples.insert_sample_id_into_db(sample_caller_connection, caller, sample_db, False)
    log.logit(f"Adding Variant IDs to the {caller} variants")
    variants.insert_variant_id_into_db(sample_caller_connection, caller, variant_db, False)
    sample_ids_block = sample_caller_connection.execute(f"SELECT DISTINCT sample_id FROM {caller}").df()["sample_id"].tolist()

    if any(sample_id in sample_ids for sample_id in sample_ids_block):
        # Create link to safe folder
        existing_sample_ids = [sample_id for sample_id in sample_ids_block if sample_id in sample_ids]
        log.logit(f"These samples: {existing_sample_ids} already exists in {caller}, need check")
        # Write a DuckDB table back to a Parquet file
        parquetPath = check_folder + os.path.basename(file).replace(".db", ".parquet")
        sample_caller_connection.execute(f"COPY {caller} TO '{parquetPath}' (FORMAT 'parquet')") 
        check = True
    else:
        # Create link to unsafe folder
        log.logit(f"These samples: {sample_ids_block} do not exist in {caller}, no check")
        # Write a DuckDB table back to a Parquet file
        parquetPath = no_check_folder + os.path.basename(file).replace(".db", ".parquet")
        sample_caller_connection.execute(f"COPY {caller} TO '{parquetPath}' (FORMAT 'parquet')")
        #sample_ids.extend(sample_ids_block)
    sample_caller_connection.close()
    return check

#def put_parquet_file_into_db_with_id_check(db_folder, db_name, parquet_folder, clobber):
def merge_caller_tables(db_path, caller_connection, variant_db, sample_db, batch_number, caller, cores, debug):
    # Setting memory_limit
    caller_connection.execute("PRAGMA memory_limit='16GB'")

    # Getting all sample names from the caller table
    sql = f"""
        SELECT DISTINCT sample_id
        FROM {caller}
    """
    caller_connection.execute(sql)
    sample_ids = caller_connection.df()["sample_id"].tolist()

    # Create no_check and check merge folders within the db_path
    check_folder = f"{db_path}/check/"
    os.makedirs(check_folder, exist_ok=True)
    check_count = 0
    no_check_folder = f"{db_path}/no_check/"
    os.makedirs(no_check_folder, exist_ok=True)
    no_check_count = 0

    # Prior to importing all parquet files, we need to insert variant_id and sample_id
    log.logit(f"Creating temporary DuckDB to insert variant_id and sample_id into parquet files")
    with indent(4, quote=' >'):
        with mp.Pool(cores) as p:
            check = p.starmap(prepare_parquet_files, [(index, db_file, caller, sample_ids, sample_db, variant_db, check_folder, no_check_folder) for index, db_file in enumerate(glob.glob(db_path + "/" + "*.db"))])
    # The multiprocessing will return either True or False depending on it needs to check or not.
    check_count = check.count(True)             # Number of files need to be checked
    no_check_count = check.count(False)         # Number of files that don't need to be checked and can perform unsafe merge            
    log.logit(f"Finished creating check and no_check folders, check: {check_count}, no_check: {no_check_count}")

    # Insert no_check folder
    if no_check_count > 0:
        log.logit(f"Inserting {no_check_count} files into {caller}")
        sql = f"INSERT INTO {caller} SELECT * FROM read_parquet('{no_check_folder}*.parquet')"
        caller_connection.execute(sql)

    for i, file in enumerate(glob.glob(check_folder + "*.parquet")):
        log.logit(f"Using safe merge for file number: {i} - {file}")
        sql = f"""
            SELECT DISTINCT sample_id
            FROM read_parquet('{file}')
        """
        sample_id = caller_connection.execute(sql).df()["sample_id"].tolist()[0]       
        sql = f"""
            INSERT INTO {caller} 
            SELECT * FROM read_parquet('{file}') 
            WHERE variant_id NOT IN (
                SELECT c.variant_id 
                FROM {caller} c 
                WHERE c.sample_id = '{sample_id}'
            )
        """
        caller_connection.execute(sql)

    total = caller_connection.execute(f"SELECT COUNT(*) FROM {caller}").fetchall()[0][0]
    caller_connection.close()
    # Remove check and no_check folders
    shutil.rmtree(check_folder)
    shutil.rmtree(no_check_folder)
    log.logit(f"Finished inserting {caller} variants, total: {total}")

def insert_caller_batch(db_path, caller_db, variant_db, sample_db, caller, batch_number, cores, debug, clobber):
    log.logit(f"Inserting variants from batch: {batch_number} into {caller_db}", color="green")
    caller_connection = db.duckdb_connect_rw(caller_db, clobber)
    setup_caller_tbl(caller_connection, caller)
    merge_caller_tables(db_path, caller_connection, variant_db, sample_db, batch_number, caller, cores, debug)
    caller_connection.close()
    log.logit(f"Finished inserting variants")
    log.logit(f"All Done!", color="green")

def annotate_fisher_test(pileup_db, caller_db, caller, batch_number, by_chromosome, debug):
    if by_chromosome:
        chromosome = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 
               'chr21', 'chr22', 'chrX', 'chrY']
    else:
        chromosome = ['ALL Chromosomes']
    log.logit(f"Performing the Fisher's Exact Test on the variants inside {caller_db} for batch: {batch_number}")
    caller = "mutect" if caller.lower() == "mutect" else "vardict"
    temp_connection = db.duckdb_connect_rw("temp_fishers.db", False)
    log.logit(f"Finding all variants within {caller_db} that does not have the fisher's exact test p-value calculated for batch: {batch_number}")
    temp_connection.execute("PRAGMA memory_limit='16GB'")
    temp_connection.execute(f"ATTACH '{pileup_db}' as pileup (READ_ONLY)")
    temp_connection.execute(f"ATTACH '{caller_db}' as caller_db (READ_ONLY)")
    for chrom in chromosome:
        log.logit(f"Processing {chrom}")
        if by_chromosome:
            filter_string = f"c.key LIKE '{chrom}:%'"
        else:
            filter_string = "TRUE"
        sql = f'''
            CREATE TABLE fisher_variants AS
            SELECT c.variant_id, c.sample_id, v.PoN_RefDepth, v.PoN_AltDepth, c.format_ref_fwd, c.format_ref_rev, c.format_alt_fwd, c.format_alt_rev
            FROM caller_db.{caller} c LEFT JOIN pileup.pileup v
            ON c.variant_id = v.variant_id
            WHERE c.fisher_p_value is NULL AND
                PoN_RefDepth is NOT NULL AND
                PoN_AltDepth is NOT NULL AND
                c.batch = {batch_number} AND
                {filter_string};
            
            SELECT * 
            FROM fisher_variants;
        '''
        if debug: log.logit(f"Executing: {sql}")
        df = temp_connection.execute(sql).df()
        temp_connection.execute(f"DROP TABLE fisher_variants;")
        if debug: log.logit(f"SQL Complete")
        length = len(df)
        log.logit(f"There were {length} variants without fisher test p-value within {caller_db}")
        if length > 0:
            log.logit(f"Calculating Fisher Exact Test for all variants inside {caller_db}")
            df = fisher_test.pvalue_df(df)
            caller_connection = db.duckdb_connect_rw(f"{caller_db}", False)
            log.logit(f"Updating {caller_db} with the fisher's exact test p-values")
            sql = f"""
                UPDATE {caller} as c
                SET fisher_p_value = df.pvalue
                FROM df
                WHERE c.variant_id = df.variant_id AND c.sample_id = df.sample_id
            """
            caller_connection.sql(sql)
            caller_connection.close()
        else:
            log.logit(f"There are no variants needed to update within {caller_db}")
    temp_connection.close()
    os.remove(f"temp_fishers.db")
    log.logit(f"Finished updating fisher test p-values inside {caller_db}")
    log.logit(f"Done!", color = "green")

def caller_to_chromosome(caller_db, caller, batch_number, chrom, base_db, debug):
    caller = "mutect" if caller.lower() == "mutect" else "vardict"
    log.logit(f"Processing {chrom}...")
    chromosome_connection = db.duckdb_connect_rw(f"{base_db}.{chrom}.db", True)
    chromosome_connection.execute("PRAGMA memory_limit='16GB'")
    chromosome_connection.execute(f"ATTACH '{caller_db}' as caller (READ_ONLY)")
    sql = f"""
            SELECT *
            FROM caller.{caller} c
            WHERE key LIKE '{chrom}:%'
        """
    if batch_number is not None:
        sql = sql + f" AND batch = {batch_number}"
    log.logit(f"Writing out variants to {base_db}.{chrom}.db")
    chromosome_connection.execute(f"CREATE TABLE {caller} AS {sql}")
    chromosome_connection.execute("DETACH caller")
    chromosome_connection.close()
    log.logit(f"Finished processing {caller_db}")
    log.logit(f"Done!", color = "green")

def merge_chromosomes(chr_path, caller_connection, caller, debug):
    # Setting memory_limit
    caller_connection.execute("PRAGMA memory_limit='16GB'")
    log.logit(f"Creating parquet files")
    with indent(4, quote=' >'):
        for i, file in enumerate(glob.glob(chr_path + "/" + "*.chr*.db")):
            log.logit(f"Processing: {i} - {file}")
            parquetPath = os.path.basename(file).replace(".db", ".parquet")
            chromosome_connection = db.duckdb_connect_rw(file, False)
            chromosome_connection.execute(f"COPY {caller} TO '{parquetPath}' (FORMAT 'parquet')")
            chromosome_connection.close()
    log.logit(f"Finished creating parquet files for all chromosomes in {chr_path}")
    log.logit(f"Inserting chromosome files into {caller}")
    sql = f"INSERT INTO {caller} SELECT * FROM read_parquet('*.parquet')"
    caller_connection.execute(sql)
    total = caller_connection.execute(f"SELECT COUNT(*) FROM {caller}").fetchall()[0][0]
    caller_connection.close()
    log.logit(f"Finished inserting {total} VCF rows into {caller}")

def chromosome_to_caller(chr_path, caller_db, caller, debug):
    log.logit(f"Inserting chromosomes from {chr_path} into {caller_db}", color="green")
    caller_connection = db.duckdb_connect_rw(caller_db, True)
    setup_caller_tbl(caller_connection, caller)
    merge_chromosomes(chr_path, caller_connection, caller, debug)
    caller_connection.close()
    log.logit(f"Finished inserting variants")
    log.logit(f"All Done!", color="green")