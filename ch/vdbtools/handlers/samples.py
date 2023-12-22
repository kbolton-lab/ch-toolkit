import os
import duckdb
import ch.vdbtools.handlers.vcf as vcf
import pandas as pd

import ch.utils.logger as log
import ch.utils.database as db

def ensure_samples_tbl(connection):
    log.logit("Ensuring or creating the samples table")
    sql = '''
        CREATE TABLE IF NOT EXISTS samples (
            sample_id          integer NOT NULL PRIMARY KEY,
            sample_name        varchar(50) NOT NULL UNIQUE,
            batch              integer NOT NULL
        )
    '''
    connection.execute(sql)

def insert_samples(samples_file, sample_duckdb, batch_number, debug, clobber):
    log.logit(f"Import Sample Information from file: {samples_file}", color="green")
    connection = db.duckdb_connect_rw(sample_duckdb, clobber)
    ensure_samples_tbl(connection)
    sample_id = connection.execute(f"SELECT MAX(sample_id) FROM samples;").fetchone()[0]
    res = pd.read_csv(samples_file,
            comment='#',
            sep='\t',
            header=None).rename(columns={0: "sample_name"})
    for sample in res['sample_name']:
        s = connection.execute(f"SELECT sample_name FROM samples WHERE sample_name = '{sample}'").fetchone()
        s = s[0] if s is not None else ""
        if sample == s:
            log.logit(f"WARNING: {sample} already exists in {sample_duckdb}.", color="yellow")
            res = res[res['sample_name'] != sample]
            res = res.reset_index(drop=True)
    sample_id = 1 if sample_id is None else sample_id + 1
    df = pd.concat([pd.Series(list(range(sample_id, len(res)+sample_id))), res], axis=1).rename(columns={0: "sample_id"})
    df['batch'] = batch_number
    vcf.duckdb_load_df_file(connection, df, "samples")
    connection.close()
    log.logit(f"Finished inserting the samples from {samples_file}")
    log.logit(f"All Done!", color="green")

def insert_sample_id_into_df(df, connection, debug):
    log.logit(f"Inserting sample_id for samples")
    connection.execute("PRAGMA memory_limit='16GB'")
    sql = f"""
            SELECT sample_id, sample_name
            FROM samples s
            WHERE s.sample_name IN (
                SELECT sample_name
                FROM df
            )
    """
    if debug: log.logit(f"Executing: {sql}")
    s = connection.execute(sql).df()
    if debug: log.logit(f"SQL Complete")
    df = s.merge(df, on='sample_name', how='left')
    df = df.drop('sample_name', axis=1)
    return df

def insert_sample_id_into_db(db, caller, sample_db, debug):
    log.logit(f"Inserting sample_id from {sample_db} for {caller} samples")
    db.execute("PRAGMA memory_limit='16GB'")
    db.execute(f"ATTACH \'{sample_db}\' as s (READ_ONLY)")
    sql = f"""
        UPDATE {caller}
        SET sample_id = s.sample_id
        FROM s.samples s
        WHERE {caller}.sample_name = s.sample_name;
    """
    if debug: log.logit(f"Executing: {sql}")
    db.execute(sql)
    db.execute(f"DETACH s")
    if debug: log.logit(f"SQL Complete")
    if debug: log.logit(f"Successfully Updated sample_id")
