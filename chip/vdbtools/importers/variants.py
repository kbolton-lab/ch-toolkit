import os, csv, glob
import vcfpy
import duckdb

import chip.vdbtools.importers.vcf as vcf
import chip.utils.logger as log
import chip.utils.database as db
import chip.utils.csv_utils as csv_utils
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
            PoN_RefDepth           INTEGER,
            PoN_AltDepth           INTEGER,
            key                    VARCHAR NOT NULL
        )
    """
    connection.execute(sql)

def drop_indexes(connection):
    log.logit("Dropping existing indexes on the variants table")
    sql = "DROP INDEX IF EXISTS chrom_pos_ref_alt_idx"
    connection.execute(sql)
    # sql = "DROP INDEX IF EXISTS variant_id_idx"
    # connection.execute(sql)
    sql = "DROP INDEX IF EXISTS batch_id_idx"
    connection.execute(sql)

def create_indexes(connection):
    log.logit("Creating new indexes on the variants table")
    with indent(4, quote=' >'):
        log.logit("Generating the chrom_pos_ref_alt_idx")
        sql = """
            CREATE UNIQUE INDEX chrom_pos_ref_alt_idx
            ON variants ( chrom, pos, ref, alt )
        """
        connection.execute(sql)
        #duckdb.sql(sql).show()

        # log.logit("Generating the variant_id_idx")
        # sql = """
        #     CREATE UNIQUE INDEX variant_id_idx
        #     ON variants ( variant_id )
        # """
        # connection.execute(sql)

        log.logit("Generating the batch_id_idx")
        sql = """
            CREATE INDEX batch_id_idx
            ON variants ( batch )
        """
        connection.execute(sql)
        #duckdb.sql(sql).show()

def setup_variants_table(connection):
    log.logit("Preparing the variant database file")
    ensure_variants_table(connection)
    drop_indexes(connection)
    return connection

def merge_variants_tables(db_path, connection, batch_number, debug):
    log.logit(f'Merging Sample Variants from {db_path}')
    with indent(4, quote=' >'):
        for i, file in enumerate(glob.glob(db_path + "/" + "*.db")):
            sample_name = os.path.basename(file).split('.')[1]
            log.logit(f"Merging: {file}")
            connection.execute(f"ATTACH \'{file}\' as sample_{i}")
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
    log.logit(f"Finished merging all tables from: {db_path}")

def ingest_variant_batch(db_path, variant_db, batch_number, debug, clobber):
    log.logit(f"Ingesting variants from batch: {batch_number} into {variant_db}", color="green")
    connection = db.duckdb_connect_rw(variant_db, clobber)
    setup_variants_table(connection)
    connection.execute("ALTER TABLE variants ADD COLUMN IF NOT EXISTS variant_id BIGINT")
    merge_variants_tables(db_path, connection, batch_number, debug)
    connection.execute("UPDATE variants SET variant_id = ROWID")
    #create_indexes(connection)
    connection.close()
    log.logit(f"Finished ingesting variants")
    log.logit(f"All Done!", color="green")

def import_sample_variants(input_vcf, variant_db, batch_number, debug, clobber):
    log.logit(f"Registering variants from file: {input_vcf}", color="green")
    connection = db.duckdb_connect_rw(variant_db, clobber)
    setup_variants_table(connection)
    counts, df = vcf.vcf_to_pd(input_vcf, "variants", batch_number, debug)
    vcf.duckdb_load_df_file(connection, df, "variants")
    connection.execute("ALTER TABLE variants ADD COLUMN IF NOT EXISTS variant_id BIGINT")
    #create_indexes(connection)
    connection.close()
    log.logit(f"Finished registering variants")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def insert_variant_id(df, connection, debug):
    log.logit(f"Inserting variant_id for all variants")
    sql = f"""
            SELECT variant_id, key
            FROM variants v
            WHERE v.key IN (
                SELECT key
                FROM df
            )
    """
    v = connection.execute(sql).df()
    df = v.merge(df, on='key', how='left')
    df = df.drop('key', axis=1)
    return df

def insert_variant_keys(df, connection, debug):
    log.logit(f"Inserting variant key for all variants")
    sql = f"""
            SELECT variant_id, key
            FROM variants v
            WHERE v.variant_id IN (
                SELECT variant_id
                FROM df
            )
    """
    keys = connection.execute(sql).df()
    df = keys.merge(df, on='variant_id', how='left')
    return df

def get_variants_from_table(connection, batch_number, chromosome):
    if chromosome != None:
        log.logit(f"Grabbing variants from batch: {batch_number} and chromosome: {chromosome} from the database")
        sql = f"SELECT variant_id, chrom, pos, ref, alt FROM variants WHERE batch = {batch_number} AND chrom = \'{chromosome}\'"
    else:
        log.logit(f"Grabbing variants from batch: {batch_number} from the database")
        sql = f"SELECT variant_id, chrom, pos, ref, alt FROM variants WHERE batch = {batch_number}"
    return connection.sql(sql)

def dump_variant_batch(variant_db, header, batch_number, chromosome, debug):
    if chromosome is None:
        log.logit(f"Dumping batch: {batch_number} variants from: {variant_db} into a VCF file", color="green")
    else:
        log.logit(f"Dumping batch: {batch_number} and chromosome: {chromosome} variants from: {variant_db} into a VCF file", color="green")
    connection = db.duckdb_connect(variant_db)
    variants = get_variants_from_table(connection, batch_number, chromosome)
    #vcf.write_variants_to_vcf(variants, header, batch_number, chromosome, debug)
    vcf.variants_to_vcf(variants, header, batch_number, chromosome, debug)
    connection.close()
    log.logit(f"Finished dumping variants into VCF file")
    log.logit(f"All Done!", color="green")

def update_pileup_variants(connection, pileup_db, debug):
    log.logit(f"Updating Pileup Information: {pileup_db} into Variants")
    with indent(4, quote=' >'):
        connection.execute(f"ATTACH \'{pileup_db}\' as pileup")
        sql = f"""
            UPDATE variants as v
            SET PoN_RefDepth = p.PoN_RefDepth, PoN_AltDepth = p.PoN_AltDepth
            FROM pileup.variants p
            WHERE v.key = p.key
        """
        connection.sql(sql)
    log.logit(f"Finished annotating variants with pileup information")

def import_pon_pileup(variant_db, pon_pileup, batch_number, debug, clobber):
    pileup_db = f"batch{batch_number}.pileup.db"
    log.logit(f"Adding pileup from batch: {batch_number} into {variant_db}", color="green")
    pileup_connection = db.duckdb_connect_rw(pileup_db, clobber)
    setup_variants_table(pileup_connection)
    counts, df = vcf.vcf_to_pd(pon_pileup, "pileup", batch_number, debug)
    vcf.duckdb_load_df_file(pileup_connection, df, "variants")
    #create_indexes(connection)
    connection.close()
    log.logit(f"Importing pileup information into {variant_db}")
    variant_connection = db.duckdb_connect_rw(variant_db, False)
    update_pileup_variants(variant_connection, pileup_db, debug)
    variant_connection.close()
    os.unlink(pileup_db)
    log.logit(f"Finished importing pileup information")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")
