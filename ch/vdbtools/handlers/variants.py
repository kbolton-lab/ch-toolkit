import os, csv, glob
import duckdb

import ch.vdbtools.handlers.vcf as vcf
import ch.utils.logger as log
import ch.utils.database as db
from clint.textui import indent, puts_err, puts

def ensure_variants_table(connection):
    log.logit("Ensuring or creating the variants table")
    sql = """
        CREATE TABLE IF NOT EXISTS variants(
            chrom                  VARCHAR(5),
            pos                    INTEGER,
            ref                    VARCHAR(255) NOT NULL,
            alt                    VARCHAR(255) NOT NULL,
            snp                    BOOLEAN,
            qc_pass                BOOLEAN,
            batch                  INTEGER,
            start                  INTEGER,
            stop                   INTEGER,
            key                    VARCHAR NOT NULL
        )
    """
    connection.execute(sql)

def ensure_pileup_table(connection):
    log.logit("Ensuring or creating the pileup table")
    sql = """
        CREATE TABLE IF NOT EXISTS pileup(
            key                    VARCHAR NOT NULL,
            PoN_RefDepth           INTEGER,
            PoN_AltDepth           INTEGER,
            batch                  INTEGER,
            variant_id             INTEGER
        )
    """
    connection.execute(sql)

def merge_variants_tables(db_path, connection, batch_number, debug):
    log.logit(f'Merging Sample Variants from {db_path} using the variant key')
    connection.execute("PRAGMA memory_limit='16GB'")
    with indent(4, quote=' >'):
        for i, file in enumerate(glob.glob(db_path + "/" + "*.db")):
            sample_name = os.path.basename(file).split('.')[1]
            log.logit(f"Merging: {i} - {file}")
            connection.execute(f"ATTACH \'{file}\' as sample_{i} (READ_ONLY)")
            sql = f"""
                INSERT INTO variants SELECT s.*
                FROM sample_{i}.variants s
                WHERE s.key NOT IN (
                    SELECT key
                    FROM variants v
                    WHERE v.key IN (
                        SELECT key
                        FROM sample_{i}.variants
                    )
                )
            """
            connection.sql(sql)
            connection.execute(f"DETACH sample_{i}")
    log.logit(f"Finished merging all tables from: {db_path}")

def insert_variant_batch(db_path, variant_db, batch_number, debug, clobber):
    log.logit(f"Inserting variants from batch: {batch_number} into {variant_db}", color="green")
    connection = db.duckdb_connect_rw(variant_db, clobber)
    ensure_variants_table(connection)
    connection.execute("ALTER TABLE variants ADD COLUMN IF NOT EXISTS variant_id BIGINT")
    merge_variants_tables(db_path, connection, batch_number, debug)
    connection.execute("UPDATE variants SET variant_id = ROWID + 1 WHERE variant_id is NULL")
    connection.close()
    log.logit(f"Finished inserting variants")
    log.logit(f"All Done!", color="green")

def import_sample_variants(input_vcf, variant_db, batch_number, debug, clobber):
    log.logit(f"Registering variants from file: {input_vcf}", color="green")
    connection = db.duckdb_connect_rw(variant_db, clobber)
    ensure_variants_table(connection)
    counts, df = vcf.vcf_to_pd(input_vcf, "variants", batch_number, debug)
    connection.execute("ALTER TABLE variants ADD COLUMN IF NOT EXISTS variant_id BIGINT")
    df = df.drop_duplicates(subset='key', keep='first')
    vcf.duckdb_load_df_file(connection, df, "variants")
    connection.close()
    log.logit(f"Finished registering variants")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def insert_variant_id_into_df(df, connection, debug):
    log.logit(f"Inserting variant_id for all variants")
    connection.execute("PRAGMA memory_limit='16GB'")
    sql = f"""
            SELECT variant_id, key
            FROM variants v
            WHERE v.key IN (
                SELECT key
                FROM df
            )
    """
    if debug: log.logit(f"Executing: {sql}")
    v = connection.execute(sql).df()
    if debug: log.logit(f"SQL Complete")
    df = v.merge(df, on='key', how='left')
    df = df.drop('key', axis=1)
    return df

def insert_variant_id_into_db(db, table, variant_db, debug):
    log.logit(f"Inserting variant_id from {variant_db} for all {table} variants")
    db.execute("PRAGMA memory_limit='16GB'")
    db.execute(f"ATTACH \'{variant_db}\' as v (READ_ONLY)")
    sql = f"""
        UPDATE {table}
        SET variant_id = v.variant_id
        FROM v.variants v
        WHERE {table}.key = v.key
        AND {table}.variant_id is NULL;
    """
    if debug: log.logit(f"Executing: {sql}")
    db.execute(sql)
    db.execute(f"DETACH v")
    if debug: log.logit(f"SQL Complete")

def insert_variant_keys(df, connection, debug):
    log.logit(f"Inserting variant key for all variants")
    connection.execute("PRAGMA memory_limit='16GB'")
    sql = f"""
            SELECT variant_id, key
            FROM variants v
            WHERE v.variant_id IN (
                SELECT variant_id
                FROM df
            )
    """
    if debug: log.logit(f"Executing: {sql}")
    keys = connection.execute(sql).df()
    if debug: log.logit(f"SQL Complete")
    df = keys.merge(df, on='variant_id', how='left')
    return df

def get_variants_from_table(connection, batch_number, chromosome):
    if batch_number != None:
        if chromosome != None:
            log.logit(f"Grabbing variants from batch: {batch_number} and chromosome: {chromosome} from the database")
            sql = f"SELECT variant_id, chrom, pos, ref, alt FROM variants WHERE batch = {batch_number} AND chrom = \'{chromosome}\'"
        else:
            log.logit(f"Grabbing variants from batch: {batch_number} from the database")
            sql = f"SELECT variant_id, chrom, pos, ref, alt FROM variants WHERE batch = {batch_number}"
    else:
        if chromosome != None:
            log.logit(f"Grabbing variants from chromosome: {chromosome} from the database")
            sql = f"SELECT variant_id, chrom, pos, ref, alt FROM variants WHERE chrom = \'{chromosome}\'"
        else:
            log.logit(f"Grabbing variants from the database")
            sql = f"SELECT variant_id, chrom, pos, ref, alt FROM variants"
    return connection.sql(sql)

def dump_variant_batch(variant_db, header, batch_number, chromosome, debug):
    if chromosome is None:
        log.logit(f"Dumping batch: {batch_number} variants from: {variant_db} into a VCF file", color="green")
    else:
        log.logit(f"Dumping batch: {batch_number} and chromosome: {chromosome} variants from: {variant_db} into a VCF file", color="green")
    connection = db.duckdb_connect_ro(variant_db)
    variants = get_variants_from_table(connection, batch_number, chromosome)
    #vcf.write_variants_to_vcf(variants, header, batch_number, chromosome, debug)
    vcf.variants_to_vcf(variants, header, batch_number, chromosome, debug)
    connection.close()
    log.logit(f"Finished dumping variants into VCF file")
    log.logit(f"All Done!", color="green")

def dump_variants_pileup(variant_db, pileup_db, header, batch_number, chromosome, debug):
    if chromosome is None:
        log.logit(f"Dumping batch: {batch_number} variants from: {variant_db} into a VCF file that needs pileup", color="green")
    else:
        log.logit(f"Dumping batch: {batch_number} and chromosome: {chromosome} variants from: {variant_db} into a VCF file that needs pileup", color="green")
    variant_connection = db.duckdb_connect_ro(variant_db)
    variant_connection.execute(f"ATTACH \'{pileup_db}\' as pileup (READ_ONLY)")
    sql = f"""
        SELECT variant_id, chrom, pos, ref, alt 
        FROM variants 
        WHERE batch = {batch_number} AND 
        variants.key NOT IN (
            SELECT key
            FROM pileup.pileup
        )
    """
    variants = variant_connection.sql(sql)
    variant_connection.execute(f"DETACH pileup")
    vcf.variants_to_vcf(variants, header, batch_number, chromosome, debug)
    variant_connection.close()
    log.logit(f"Finished dumping variants that need pileup into VCF file")
    log.logit(f"All Done!", color="green")

def import_pon_pileup(pileup_db, variant_db, pon_pileup, batch_number, debug, clobber):
    log.logit(f"Adding pileup from batch: {batch_number} into {pileup_db}", color="green")
    pileup_connection = db.duckdb_connect_rw(pileup_db, clobber)
    ensure_pileup_table(pileup_connection)
    counts, df = vcf.vcf_to_pd(pon_pileup, "pileup", batch_number, debug)
    vcf.duckdb_load_df_file(pileup_connection, df, "pileup")
    log.logit(f"Adding Variant IDs to the pileup database")
    insert_variant_id_into_db(pileup_connection, "pileup", variant_db, debug)
    pileup_connection.close()
    log.logit(f"Finished importing pileup information")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def variant_to_chromosome(variant_db, variant, chrom, base_db, debug):
    log.logit(f"Processing {chrom}...")
    chromosome_connection = db.duckdb_connect_rw(f"{base_db}.{chrom}.db", True)
    chromosome_connection.execute("PRAGMA memory_limit='16GB'")
    chromosome_connection.execute(f"ATTACH '{variant_db}' as {variant} (READ_ONLY)")
    sql = f"""
            SELECT *
            FROM {variant}.variants v
            WHERE key LIKE '{chrom}:%'
        """
    log.logit(f"Writing out variants to {base_db}.{chrom}.db")
    chromosome_connection.execute(f"CREATE TABLE variants AS {sql}")
    chromosome_connection.execute(f"DETACH {variant}")
    chromosome_connection.close()
    log.logit(f"Finished processing {variant_db}")
    log.logit(f"Done!", color = "green")

def pileup_to_chromosome(pileup_db, pileup, chrom, base_db, debug):
    log.logit(f"Processing {chrom}...")
    chromosome_connection = db.duckdb_connect_rw(f"{base_db}.{chrom}.db", True)
    chromosome_connection.execute("PRAGMA memory_limit='16GB'")
    chromosome_connection.execute(f"ATTACH '{pileup_db}' as {pileup} (READ_ONLY)")
    sql = f"""
            SELECT *
            FROM {pileup}.{pileup} p
            WHERE key LIKE '{chrom}:%'
        """
    log.logit(f"Writing out variants to {base_db}.{chrom}.db")
    chromosome_connection.execute(f"CREATE TABLE {pileup} AS {sql}")
    chromosome_connection.execute(f"DETACH {pileup}")
    chromosome_connection.close()
    log.logit(f"Finished processing {pileup_db}")
    log.logit(f"Done!", color = "green")