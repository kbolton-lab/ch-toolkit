import chip.utils.logger as log
import chip.utils.redis as redis

import duckdb
from clint.textui import indent, puts_err, puts

def _ensure_variants_table(connection):
    log.logit("Ensuring or creating the variants table")
    sql = """
        CREATE TABLE IF NOT EXISTS variants(
			variant_id             BIGINT PRIMARY KEY,
			chrom                  VARCHAR(5),
			pos                    INTEGER,
			ref                    VARCHAR(255) NOT NULL,
			alt                    VARCHAR(255) NOT NULL,
			qc_pass                BOOLEAN,
			batch_id               INTEGER,
            start                  INTEGER,
            stop                   INTEGER
        )
    """
    connection.execute(sql)

def _drop_indexes(connection):
    log.logit("Dropping existing indexes on the variants table")
    sql = "DROP INDEX IF EXISTS chrom_pos_ref_alt_idx"
    connection.execute(sql)
    sql = "DROP INDEX IF EXISTS variant_id_idx"
    connection.execute(sql)
    sql = "DROP INDEX IF EXISTS batch_id_idx"
    connection.execute(sql)

def _create_indexes(connection):
    log.logit("Creating new indexes on the variants table")
    with indent(4, quote=' >'):
        log.logit("Generating the chrom_pos_ref_alt_idx")
        sql = """
            CREATE UNIQUE INDEX chrom_pos_ref_alt_idx
            ON variants ( chrom, pos, ref, alt )
        """
        connection.execute(sql)

        log.logit("Generating the variant_id_idx")
        sql = """
            CREATE UNIQUE INDEX variant_id_idx
            ON variants ( variant_id )
        """
        connection.execute(sql)

        log.logit("Generating the batch_id_idx")
        sql = """
            CREATE UNIQUE INDEX batch_id_idx
            ON variants ( batch_id )
        """
        connection.execute(sql)

def _setup_variants_table(connection):
    log.logit("Preparing the variant database file")
    _ensure_variants_table(connection)
    _drop_indexes(connection)
    return connection

def _bulk_insert(duckdb_connection, entries):
    sql = "insert into variants ( variant_id, chrom, pos, ref, alt, batch_id, start, stop ) values ( ?, ?, ?, ?, ?, ?, ?, ?)"
    log.logit(f"Starting to bulk insert variant batch into duckdb ( {len(entries)} items)")
    duckdb.executemany(sql, entries, duckdb_connection)
    log.logit(f"Finished DuckDB insertion")

def _insert_variants(redis_db, duckdb_connection, batch_number, debug):
    log.logit(f"Starting to ingest variants in batch: {batch_number}")
    redis_set_key = f"batch:{batch_number}"
    total_variants = int(redis_db.scard(redis_set_key))
    log.logit(f"Total # of variants to insert in batch {batch_number}: {total_variants}")
    log.logit("Inserting batch variants into duckdb variants table")
    variant_batch = []
    with indent(4, quote=' >'):
        for (i, key) in enumerate(redis_db.sscan_iter(redis_set_key)):
            variant_id = int(redis_db.get(key))
            (chrom, pos, ref, alt) = key.decode().split(':')
            pos = int(pos)
            start = pos
            stop = start + len(alt)
            item = (variant_id, chrom, pos, ref, alt, batch_number, start, stop)
            variant_batch.append( item )
            if debug: log.logit(f"{i} -- {item}")
            if i % 25000 == 0 and i != 0:
#            if i % 1 == 0 and i != 0:
#                import pdb; pdb.set_trace()
                _bulk_insert(duckdb_connection, variant_batch)
                log.logit(f"# Variants Processed: {i}")
                variant_batch = []

        if len(variant_batch) > 0:
            _bulk_insert(duckdb_connection, variant_batch)

    return { 'counts': total_variants }

def ingest_variant_batch(duckdb_file, redis_host, redis_port, batch_number, debug):
    log.logit(f"Ingesting variants from redis batch: {batch_number} into {duckdb_file}", color="green")
    redis_db = redis.redis_connect(redis_host, redis_port)
    duckdb_connection = duckdb.connect(duckdb_file)
    _setup_variants_table(duckdb_connection)
    counts = _insert_variants(redis_db, duckdb_connection, batch_number, debug)
    _create_indexes(duckdb_connection)
    duckdb_connection.close()
    log.logit(f"Finished ingesting variants")
    log.logit(f"Variants Processed - Total: {counts['total']}", color="green")
    log.logit(f"All Done!", color="green")
