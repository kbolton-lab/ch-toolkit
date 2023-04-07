import os
import vcfpy
import duckdb

import chip.vdbtools.importers.vcf as vcf
import chip.utils.logger as log
import chip.utils.database as db
import chip.utils.fisher_exact_test as fisher_test

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
            fisher_p_value                      decimal(22,20),
            PRIMARY KEY( sample_id, variant_id )
        )
    '''
    connection.execute(sql)

def ensure_vardict_tbl(connection):
    log.logit("Ensuring or creating the vardict table")

    sql = '''
        CREATE TYPE IF NOT EXISTS vardict_filter_type AS ENUM (
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
            fisher_p_value          decimal(22,20),
            PRIMARY KEY( sample_id, variant_id )
        )
    '''
    connection.execute(sql)

def drop_indexes_mutect(connection):
    log.logit("Dropping existing indexes on the mutect table")
    sql = "DROP INDEX IF EXISTS sample_variant_idx"
    connection.execute(sql)

def drop_indexes_vardict(connection):
    log.logit("Dropping existing indexes on the vardict table")
    sql = "DROP INDEX IF EXISTS sample_variant_idx"
    connection.execute(sql)

def create_indexes_mutect(connection):
    log.logit("Creating new indexes on the mutect table")
    log.logit("Generating the sample_variant_idx")
    sql = """
        CREATE UNIQUE INDEX sample_variant_idx
        ON mutect ( sample_id, variant_id )
    """
    connection.execute(sql)

def create_indexes_vardict(connection):
    log.logit("Creating new indexes on the vardict table")
    log.logit("Generating the sample_variant_idx")
    sql = """
        CREATE UNIQUE INDEX sample_variant_idx
        ON vardict ( sample_id, variant_id )
    """
    connection.execute(sql)

def setup_mutect_tbl(connection):
    log.logit("Preparing the mutect database file")
    ensure_mutect_tbl(connection)
    drop_indexes_mutect(connection)
    return connection

def setup_vardict_tbl(connection):
    log.logit("Preparing the vardict database file")
    ensure_vardict_tbl(connection)
    drop_indexes_vardict(connection)
    return connection

# def bulk_insert_mutect(connection, entries):
#     #sql = "INSERT INTO mutect (sample_id, variant_id, version, mutect_filter, info_as_filterstatus, info_as_sb_table, info_dp, info_ecnt, info_mbq_ref, info_mbq_alt, info_mfrl_ref, info_mfrl_alt, info_mmq_ref, info_mmq_alt, info_mpos, info_popaf, info_roq, info_rpa_ref, info_rpa_alt, info_ru, info_str, info_strq, info_tlod, format_gt, format_af, format_dp, format_ref_count, format_alt_count, format_ref_f1r2, format_alt_f1r2, format_ref_f2r1, format_alt_f2r1, format_ref_fwd, format_ref_rev, format_alt_fwd, format_alt_rev, fisher_p_value) VALUES (?, ?, '?', ?, ?, '?', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, '?', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
#     sql = "INSERT INTO mutect (sample_id, variant_id, version, mutect_filter, info_as_filterstatus, info_as_sb_table, info_dp, info_ecnt, info_mbq_ref, info_mbq_alt, info_mfrl_ref, info_mfrl_alt, info_mmq_ref, info_mmq_alt, info_mpos, info_popaf, info_roq, info_rpa_ref, info_rpa_alt, info_ru, info_str, info_strq, info_tlod, format_gt, format_af, format_dp, format_ref_count, format_alt_count, format_ref_f1r2, format_alt_f1r2, format_ref_f2r1, format_alt_f2r1, format_ref_fwd, format_ref_rev, format_alt_fwd, format_alt_rev, fisher_p_value) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
#     log.logit(f"Starting to bulk insert mutect batch into duckdb ( {len(entries)} items)")
#     duckdb.executemany(sql, entries, connection)
#     log.logit(f"Finished DuckDB insertion")
#
# def bulk_insert_vardict(connection, entries):
#     sql = "INSERT INTO vardict (sample_id, variant_id, version, vardict_filter, info_type, info_vd, info_af, info_bias, info_refbias, info_varbias, info_pmean, info_pstd, info_qual, info_qstd, info_sbf, info_oddratio, info_mq, info_sn, info_hiaf, info_adjaf, info_shift3, info_msi, info_msilen, info_nm, info_lseq, info_rseq, info_hicnt, info_hicov, info_splitread, info_spanpair, info_duprate, format_gt, format_af, format_dp, format_ref_count, format_alt_count, format_vd, format_ref_fwd, format_ref_rev, format_alt_fwd, format_alt_rev, fisher_p_value) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
#     log.logit(f"Starting to bulk insert vardict batch into duckdb ( {len(entries)} items)")
#     duckdb.executemany(sql, entries, connection)
#     log.logit(f"Finished DuckDB insertion")
#
# def bulk_update_caller_pvalue(connection, entries, caller):
#     sql = f"UPDATE {caller} SET fisher_p_value=? WHERE variant_id=? AND sample_id=?"
#     log.logit(f"Starting to bulk update p-values into duckdb ( {len(entries)} items)")
#     duckdb.executemany(sql, entries, connection)
#     log.logit(f"Finished DuckDB update")

def process_mutect(df, debug):
    #df[['info_as_filterstatus','info_as_sb_table','info_dp','info_ecnt','info_mbq','info_mfrl','info_mmq','info_mpos','info_popaf','info_roq','info_rpa','info_ru','info_str','info_strq','info_tlod']] = df['info'].str.split(";", expand=True)
    test = df['info'].str.extract('AS_FilterStatus=(?P<info_as_filterstatus>.+?);AS_SB_TABLE=(?P<info_as_sb_table>.+?);DP=(?P<info_dp>.+?);ECNT=(?P<info_ecnt>.+?);MBQ=(?P<info_mbq>.+?)')
    print(test)
    #df = df[['sample_id', 'variant_id', 'version', 'mutect_filter_type', 'info_as_filterstatus', 'info_as_sb_table', 'info_dp', 'info_ecnt', 'info_mbq_ref', 'info_mbq_alt', 'info_mfrl_ref', 'info_mfrl_alt', 'info_mmq_ref', 'info_mmq_alt', 'info_mpos', 'info_popaf', 'info_roq', 'info_rpa_ref', 'info_rpa_alt', 'info_ru', 'info_str', 'info_strq', 'info_tlod', 'format_gt', 'format_af', 'format_dp', 'format_ref_count', 'format_alt_count', 'format_ref_f1r2', 'format_alt_f1r2', 'format_ref_f2r1', 'format_alt_f2r1', 'format_ref_fwd', 'format_ref_rev', 'format_alt_fwd', 'format_alt_rev', 'fisher_p_value']]
    return df

def insert_mutect_caller_(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug):
    log.logit(f"Registering mutect variants from: {input_vcf} in batch: {batch_number}", color="green")
    sample_name = os.path.basename(input_vcf).split('.')[1]
    log.logit(f"Grabbing sample_id from {sample_db} for sample: {sample_name}")
    sample_id = db.duckdb_connect(sample_db).execute(f"SELECT sample_id FROM samples WHERE sample_name = '{sample_name}';").fetchone()[0]
    caller_connection = db.duckdb_connect_rw(db_path, clobber)
    setup_mutect_tbl(caller_connection)
    counts, df = vcf.vcf_to_pd(input_vcf, "caller", batch_number, debug)
    df = process_mutect(df, debug)
    print(df)
    #variant_connection = db.duckdb_connect_ro(variant_db)
    #insert_variant_id(df, variant_connection, debug)
    #vcf.duckdb_load_df_file(connection, df, "mutect")
    #create_indexes(connection)
    #connection.close()

# Vardict
# res = res [['sample_id', 'variant_id', 'version', 'filter_type', 'info_type', 'info_vd', 'info_af', 'info_bias', 'info_refbias', 'info_varbias', 'info_pmean', 'info_pstd', 'info_qual', 'info_qstd', 'info_sbf', 'info_oddratio', 'info_mq', 'info_sn', 'info_hiaf', 'info_adjaf', 'info_shift3', 'info_msi', 'info_msilen', 'info_nm', 'info_lseq', 'info_rseq', 'info_hicnt', 'info_hicov', 'info_splitread', 'info_spanpair', 'info_duprate', 'format_gt', 'format_af', 'format_dp', 'format_ref_count', 'format_alt_count', 'format_vd', 'format_ref_fwd', 'format_ref_rev', 'format_alt_fwd', 'format_alt_rev', 'fisher_p_value']]
# filter_type = record.FILTER

# info_type = record.INFO.get('TYPE')
# info_dp = record.INFO.get('DP', 0)
# info_vd = record.INFO.get('VD', 0)
# info_af = record.INFO.get('AF', 0)[0]
# info_bias = record.INFO.get('BIAS')
# info_refbias = record.INFO.get('REFBIAS')
# info_varbias = record.INFO.get('VARBIAS')
# info_pmean = record.INFO.get('PMEAN')
# info_pstd = record.INFO.get('PSTD')
# info_qual = record.INFO.get('QUAL')
# info_qstd = record.INFO.get('QSTD')
# info_sbf = record.INFO.get('SBF')
# info_oddratio = record.INFO.get('ODDRATIO')
# info_mq = record.INFO.get('MQ')
# info_sn = record.INFO.get('SN')
# info_hiaf = record.INFO.get('HIAF')
# info_adjaf = record.INFO.get('ADJAF')
# info_shift3 = record.INFO.get('SHIFT3')
# info_msi = record.INFO.get('MSI')
# info_msilen = record.INFO.get('MSILEN')
# info_nm = record.INFO.get('NM')
# info_lseq = record.INFO.get('LSEQ')
# info_rseq = record.INFO.get('RSEQ')
# info_hicnt = record.INFO.get('HICNT')
# info_hicov = record.INFO.get('HICOV')
# info_splitread = record.INFO.get('SPLITREAD')
# info_spanpair = record.INFO.get('SPANPAIR')
# info_duprate = record.INFO.get('DUPRATE')
#
# format_gt = record.calls[0].data.get('GT')
# format_dp = record.calls[0].data.get('DP', 0)
# format_af = record.calls[0].data.get('AF')[0]
# (format_ref_count, format_alt_count) = record.calls[0].data.get('AD')
# (format_ref_fwd, format_ref_rev) = record.calls[0].data.get('RD')
# (format_alt_fwd, format_alt_rev) = record.calls[0].data.get('ALD')
# format_vd = record.calls[0].data.get('VD', 0)
# fisher_p_value = fisher_test.pvalue(PoN_RefDepth, PoN_AltDepth, format_ref_fwd+format_ref_rev, format_alt_fwd+format_alt_rev)

def insert_mutect_caller(db_path, input_vcf, chromosome, variant_duckdb, sample_duckdb, clobber, debug):
    con = db.duckdb_connect_rw(db_path, clobber)
    setup_mutect_tbl(con)
    reader = vcfpy.Reader.from_path(input_vcf)
    #version = reader.header.HeaderLine("MutectVersion")    # No Idea how to get this information yet...
    version="2.2"
    window_size = 10_000
    sample_name = os.path.basename(input_vcf).split('.')[1]
    variants = reader.fetch(chromosome) if chromosome != None else reader
    variant_duckdb_connection = db.duckdb_connect_ro(variant_duckdb)
    sample_duckdb_connection = db.duckdb_connect(sample_duckdb)
    window = []
    count = 1
    # Get the Sample_ID
    sample_id = sample_duckdb_connection.execute(f"SELECT sample_id FROM samples WHERE sample_name = '{sample_name}';").fetchone()[0]
    for record in variants:
        # Get the Variant_ID
        variant_id, PoN_RefDepth, PoN_AltDepth = variant_duckdb_connection.execute(f"SELECT variant_id, PoN_RefDepth, PoN_AltDepth FROM variants WHERE chrom = '{record.CHROM}' AND pos = {record.POS} AND ref = '{record.REF}' AND alt = '{record.ALT[0].value}';").fetchone()
        if debug: log.logit(f"{variant_id} - {PoN_RefDepth}, {PoN_AltDepth}")

        mutect_filter_type = record.FILTER
        info_as_filterstatus = ['Multiallelic'] if record.INFO.get('AS_FilterStatus') is None else record.INFO.get('AS_FilterStatus')
        info_as_sb_table = record.INFO.get('AS_SB_TABLE')
        info_dp = record.INFO.get('DP', 0)
        info_ecnt = record.INFO.get('ECNT', 0)
        info_mbq_ref = record.INFO.get("MBQ")[0]
        info_mbq_alt = record.INFO.get("MBQ")[1]
        info_mfrl_ref = record.INFO.get("MFRL")[0]
        info_mfrl_alt = record.INFO.get("MFRL")[1]
        info_mmq_ref = record.INFO.get("MMQ")[0]
        info_mmq_alt = record.INFO.get("MMQ")[1]
        info_mpos = record.INFO.get("MPOS")[0]
        info_popaf = record.INFO.get("POPAF")[0]
        info_roq = record.INFO.get("ROQ")
        info_rpa_ref = None if record.INFO.get("RPA") is None else record.INFO.get("RPA")[0]
        info_rpa_alt = None if record.INFO.get("RPA") is None else record.INFO.get("RPA")[1]
        info_ru = None if record.INFO.get("RU") is None else record.INFO.get("RU")[0]
        info_str = record.INFO.get("STR", False)
        info_strq = record.INFO.get("STRQ", None)
        info_tlod = record.INFO.get("TLOD")[0]

        format_dp = record.calls[0].data.get('DP', 0)
        format_af = record.calls[0].data.get('AF')[0]
        (format_ref_count, format_alt_count) = record.calls[0].data.get('AD')
        (format_ref_f1r2, format_alt_f1r2) = record.calls[0].data.get('F1R2')
        (format_ref_f2r1, format_alt_f2r1) = record.calls[0].data.get('F2R1')
        format_gt = record.calls[0].data.get('GT')
        (format_ref_fwd, format_ref_rev, format_alt_fwd, format_alt_rev) = record.calls[0].data.get('SB')
        fisher_p_value = fisher_test.pvalue(PoN_RefDepth, PoN_AltDepth, format_ref_fwd+format_ref_rev, format_alt_fwd+format_alt_rev)

        if debug: log.logit(f"sample_id={sample_id} | variant_id={variant_id} | version={version} | mutect_filter_type={mutect_filter_type} \nINFO: info_as_filterstatus={info_as_filterstatus} | info_as_sb_table={info_as_sb_table} | info_dp={info_dp} | info_ecnt={info_ecnt} | info_mbq_ref={info_mbq_ref} | info_mbq_alt={info_mbq_alt} | info_mfrl_ref={info_mfrl_ref} | info_mfrl_alt={info_mfrl_alt} | info_mmq_ref={info_mmq_ref} | info_mmq_alt={info_mmq_alt} | info_mpos={info_mpos} | info_popaf={info_popaf} | info_roq={info_roq} | info_rpa_ref={info_rpa_ref} | info_rpa_alt={info_rpa_alt} | info_ru={info_ru} | info_str={info_str} | info_strq={info_strq} | info_tlod={info_tlod} \nFORMAT: format_dp={format_dp} | format_af={format_af} | format_ref_count={format_ref_count} | format_alt_count={format_alt_count} | format_ref_f1r2={format_ref_f1r2} | format_alt_f1r2={format_alt_f1r2} | format_ref_f2r1={format_ref_f2r1} | format_alt_f2r1={format_alt_f2r1} | format_gt={format_gt} | format_ref_fwd={format_ref_fwd} | format_ref_rev={format_ref_rev} | format_alt_fwd={format_alt_fwd} | format_alt_rev={format_alt_rev} | fisher_p_value={fisher_p_value}")

        row = (sample_id, variant_id, version, mutect_filter_type, info_as_filterstatus, info_as_sb_table, info_dp, info_ecnt, info_mbq_ref, info_mbq_alt, info_mfrl_ref, info_mfrl_alt, info_mmq_ref, info_mmq_alt, info_mpos, info_popaf, info_roq, info_rpa_ref, info_rpa_alt, info_ru, info_str, info_strq, info_tlod, format_gt, format_af, format_dp, format_ref_count, format_alt_count, format_ref_f1r2, format_alt_f1r2, format_ref_f2r1, format_alt_f2r1, format_ref_fwd, format_ref_rev, format_alt_fwd, format_alt_rev, fisher_p_value)
        window.append(row)

        if count % window_size == 0:
            bulk_insert_mutect(con, window)
            log.logit(f"# Variants Processed: {count}", color = "green")
            window = []
        count += 1
    if len(window) > 0:
        bulk_insert_mutect(con, window)
    variant_duckdb_connection.close()
    create_indexes_mutect(con)
    con.close()
    log.logit(f"Finished inserting mutect variants")
    log.logit(f"Variants Processed - Total: {count}", color="green")
    log.logit(f"All Done!", color="green")

def insert_vardict_caller(db_path, input_vcf, chromosome, variant_duckdb, sample_duckdb, clobber, debug):
    con = db.duckdb_connect_rw(db_path, clobber)
    setup_vardict_tbl(con)
    reader = vcfpy.Reader.from_path(input_vcf)
    version="v1.8.2"

    sample_name = os.path.basename(input_vcf).split('.')[1]
    variants = reader.fetch(chromosome) if chromosome != None else reader
    variant_duckdb_connection = db.duckdb_connect_ro(variant_duckdb)
    sample_duckdb_connection = db.duckdb_connect(sample_duckdb)
    window = []
    count = 1
    for record in variants:
        # Get the Variant_ID
        variant_id, PoN_RefDepth, PoN_AltDepth = variant_duckdb_connection.execute(f"SELECT variant_id, PoN_RefDepth, PoN_AltDepth FROM variants WHERE chrom = '{record.CHROM}' AND pos = {record.POS} AND ref = '{record.REF}' AND alt = '{record.ALT[0].value}';").fetchone()
        if debug: log.logit(f"{variant_id} - {PoN_RefDepth}, {PoN_AltDepth}")
        # Get the Sample_ID
        sample_id = sample_duckdb_connection.execute(f"SELECT sample_id FROM samples WHERE sample_name = '{sample_name}';").fetchone()[0]

        filter_type = record.FILTER
        info_type = record.INFO.get('TYPE')
        info_dp = record.INFO.get('DP', 0)
        info_vd = record.INFO.get('VD', 0)
        info_af = record.INFO.get('AF', 0)[0]
        info_bias = record.INFO.get('BIAS')
        info_refbias = record.INFO.get('REFBIAS')
        info_varbias = record.INFO.get('VARBIAS')
        info_pmean = record.INFO.get('PMEAN')
        info_pstd = record.INFO.get('PSTD')
        info_qual = record.INFO.get('QUAL')
        info_qstd = record.INFO.get('QSTD')
        info_sbf = record.INFO.get('SBF')
        info_oddratio = record.INFO.get('ODDRATIO')
        info_mq = record.INFO.get('MQ')
        info_sn = record.INFO.get('SN')
        info_hiaf = record.INFO.get('HIAF')
        info_adjaf = record.INFO.get('ADJAF')
        info_shift3 = record.INFO.get('SHIFT3')
        info_msi = record.INFO.get('MSI')
        info_msilen = record.INFO.get('MSILEN')
        info_nm = record.INFO.get('NM')
        info_lseq = record.INFO.get('LSEQ')
        info_rseq = record.INFO.get('RSEQ')
        info_hicnt = record.INFO.get('HICNT')
        info_hicov = record.INFO.get('HICOV')
        info_splitread = record.INFO.get('SPLITREAD')
        info_spanpair = record.INFO.get('SPANPAIR')
        info_duprate = record.INFO.get('DUPRATE')

        format_gt = record.calls[0].data.get('GT')
        format_dp = record.calls[0].data.get('DP', 0)
        format_af = record.calls[0].data.get('AF')[0]
        (format_ref_count, format_alt_count) = record.calls[0].data.get('AD')
        (format_ref_fwd, format_ref_rev) = record.calls[0].data.get('RD')
        (format_alt_fwd, format_alt_rev) = record.calls[0].data.get('ALD')
        format_vd = record.calls[0].data.get('VD', 0)
        fisher_p_value = fisher_test.pvalue(PoN_RefDepth, PoN_AltDepth, format_ref_fwd+format_ref_rev, format_alt_fwd+format_alt_rev)

        if debug: log.logit(f"sample_id={sample_id} | variant_id={variant_id} | version={version} | vardict_filter_type={filter_type} \nINFO: info_type={info_type} | info_vd={info_vd} | info_af={info_af} | info_bias={info_bias} | info_refbias={info_refbias} | info_varbias={info_varbias} | info_pmean={info_pmean} | info_pstd={info_pstd} | info_qual={info_qual} | info_qstd={info_qstd} | info_sbf={info_sbf} | info_oddratio={info_oddratio} | info_mq={info_mq} | info_sn={info_sn} | info_hiaf={info_hiaf} | info_adjaf={info_adjaf} | info_shift3={info_shift3} | info_msi={info_msi} | info_msilen={info_msilen} | info_nm={info_nm} | info_lseq={info_lseq} | info_rseq={info_rseq} | info_hicnt={info_hicnt} | info_hicov={info_hicov} | info_splitread={info_splitread} | info_spanpair={info_spanpair} | info_duprate={info_duprate} \nFORMAT: format_dp={format_dp} | format_af={format_af} | format_ref_count={format_ref_count} | format_alt_count={format_alt_count} | format_vd={format_vd} | format_gt={format_gt} | format_ref_fwd={format_ref_fwd} | format_ref_rev={format_ref_rev} | format_alt_fwd={format_alt_fwd} | format_alt_rev={format_alt_rev} | fisher_p_value={fisher_p_value}")

        row = (sample_id, variant_id, version, filter_type, info_type, info_vd, info_af, info_bias, info_refbias, info_varbias, info_pmean, info_pstd, info_qual, info_qstd, info_sbf, info_oddratio, info_mq, info_sn, info_hiaf, info_adjaf, info_shift3, info_msi, info_msilen, info_nm, info_lseq, info_rseq, info_hicnt, info_hicov, info_splitread, info_spanpair, info_duprate, format_gt, format_af, format_dp, format_ref_count, format_alt_count, format_vd, format_ref_fwd, format_ref_rev, format_alt_fwd, format_alt_rev, fisher_p_value)
        window.append(row)

        if count % window_size == 0:
            bulk_insert_vardict(con, window)
            log.logit(f"# Variants Processed: {count}", color = "green")
            window = []
        count += 1
    if len(window) > 0:
        bulk_insert_vardict(con, window)
    variant_duckdb_connection.close()
    create_indexes_vardict(con)
    con.close()
    log.logit(f"Finished inserting vardict variants")
    log.logit(f"Variants Processed - Total: {count}", color="green")
    log.logit(f"All Done!", color="green")

def annotate_fisher_test(variant_duckdb, caller_duckdb, caller, batch_number, chromosome, window_size, debug):
    variant_duckdb_connection = db.duckdb_connect_rw(variant_duckdb, False)
    caller_duckdb_connection = db.duckdb_connect_rw(caller_duckdb, False)

    caller = "mutect" if caller.lower() == "mutect" else "vardict"

    setup_mutect_tbl(variant_duckdb_connection)
    create_indexes_mutect(variant_duckdb_connection)
    log.logit(f"Finding all variants within {caller_duckdb} that does not have the fisher pvalue calcualted...")
    variant_duckdb_connection.execute(f"ATTACH '{caller_duckdb}' as {caller}")
    variant_duckdb_connection.execute(f"INSERT INTO {caller} SELECT * FROM {caller}.{caller} WHERE {caller}.fisher_p_value is NULL")
    sql = f'''
        SELECT caller.variant_id, caller.sample_id, caller.format_ref_fwd, caller.format_ref_rev, caller.format_alt_fwd, caller.format_alt_rev, variants.PoN_RefDepth, variants.PoN_AltDepth
        FROM {caller} caller LEFT JOIN variants
        ON caller.variant_id = variants.variant_id
    '''
    res = variant_duckdb_connection.execute(sql)
    from itertools import starmap
    log.logit(f"Calculating Fisher Exact Test for all variants inside {caller_duckdb}")
    res = list(starmap(fisher_test.pvalue_variant_sample, res.fetchall()))
    length = len(res)
    log.logit(f"There were {length} variants without fisher test p-value within {caller_duckdb}")
    log.logit(f"Removing the attached tables in {variant_duckdb}")
    variant_duckdb_connection.execute(f"DROP TABLE {caller}")
    variant_duckdb_connection.execute(f"DROP type {caller}_filter_type")
    variant_duckdb_connection.close()
    window = []
    log.logit(f"Annotating {caller_duckdb} with the p-values")
    for i, row in enumerate(res):
        # pvalue = row[0]
        # id = row[1]
        # sample = row[2]
        # if debug: log.logit(f"{i} -- {pvalue}, {id}, {sample}")
        # sql = f"UPDATE {caller} SET fisher_p_value={pvalue} WHERE variant_id={id} AND sample_id={sample}"
        # if debug: log.logit(f"{sql}")
        # caller_duckdb_connection.execute(sql)
        window.append(row)
        if debug: log.logit(f"{i} -- {row}")
        if i % window_size == 0 and i != 0:
            bulk_update_caller_pvalue(caller_duckdb_connection, window, caller)
            log.logit(f"# Variants Processed: {i}")
            window = []
    if len(window) > 0:
        bulk_update_caller_pvalue(caller_duckdb_connection, window, caller)
    caller_duckdb_connection.close()
    log.logit(f"Finished updating fisher test pvalues inside {caller_duckdb}")
    log.logit(f"Done!", color = "green")
