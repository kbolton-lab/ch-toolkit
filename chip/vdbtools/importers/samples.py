import os
import duckdb
import chip.vdbtools.importers.vcf as vcf
import pandas as pd

import chip.utils.logger as log
import chip.utils.database as db

def ensure_samples_tbl(connection):
    log.logit("Ensuring or creating the mutect table")
    sql = '''
        CREATE TABLE IF NOT EXISTS samples (
            sample_id          integer NOT NULL PRIMARY KEY,
            sample_name        varchar(50) NOT NULL UNIQUE
        )
    '''
    connection.execute(sql)

def drop_indexes_samples(connection):
    log.logit("Dropping existing indexes on the samples table")
    sql = "DROP INDEX IF EXISTS sample_id_idx"
    connection.execute(sql)

def create_indexes_samples(connection):
    log.logit("Creating new indexes on the samples table")
    log.logit("Generating the sample_id_idx")
    sql = """
        CREATE UNIQUE INDEX sample_id_idx
        ON samples ( sample_id )
    """
    connection.execute(sql)

def setup_samples_tbl(connection):
    log.logit("Preparing the samples database file")
    ensure_samples_tbl(connection)
    drop_indexes_samples(connection)
    return connection

def insert_samples(samples_file, sample_duckdb, debug, clobber):
    connection = db.duckdb_connect_rw(sample_duckdb, clobber)
    setup_samples_tbl(connection)
    sample_id = connection.execute(f"SELECT MAX(sample_id) FROM samples;").fetchone()[0]
    sample_id = 1 if sample_id is None else sample_id + 1
    res = pd.read_csv(samples_file,
            comment='#',
            sep='\t',
            header=None).rename(columns={0: "sample_name"})
    df = pd.concat([pd.Series(list(range(sample_id, len(res)+sample_id))), res], axis=1).rename(columns={0: "sample_id"})
    vcf.duckdb_load_df_file(connection, df, "samples")
    create_indexes_samples(connection)
    connection.close()
    log.logit(f"Finished inserting the samples from {samples_file}")
    log.logit(f"All Done!", color="green")
