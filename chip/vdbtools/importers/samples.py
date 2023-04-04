import csv
import os
import duckdb

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

def load_sample_csv_file(connection, csv):
    sql = f"""
        COPY variants FROM '{csv}' (AUTO_DETECT TRUE)
    """
    log.logit(f"Starting to load csv into duckdb")
    duckdb_connection.execute(sql)
    log.logit(f"Finished loading csv into duckdb")

def bulk_insert_samples(connection, entries):
    sql = f"INSERT INTO samples (sample_id, sample_name) VALUES (?, ?)"
    log.logit(f"Starting to bulk insert samples into duckdb ( {len(entries)} items)")
    duckdb.executemany(sql, entries, connection)
    log.logit(f"Finished Samples Insertion")

def insert_samples(samples, sample_duckdb, debug, clobber):
    con = db.duckdb_connect_rw(sample_duckdb, clobber)
    setup_samples_tbl(con)
    sample_id = con.execute(f"SELECT MAX(sample_id) FROM samples;").fetchone()[0]
    sample_id = 1 if sample_id is None else sample_id + 1
    samples_reader = open(samples, 'r')
    window = []
    for line in samples_reader.readlines():
        row = (sample_id , line.strip("\n"))
        if debug: log.logit(f"{row}")
        window.append(row)
        sample_id += 1
    samples_reader.close()
    bulk_insert_samples(con, window)
    create_indexes_samples(con)
    con.close()
    log.logit(f"Finished inserting the samples from {samples}")
    log.logit(f"All Done!", color="green")
