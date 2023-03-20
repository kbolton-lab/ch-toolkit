import os, gzip, csv

import chip.utils.logger as log
import chip.utils.database as db

import duckdb
from clint.textui import indent, puts_err, puts

def ensure_variants_table(connection):
    log.logit("Ensuring or creating the variants table")
    sql = """
        CREATE TABLE IF NOT EXISTS variants(
            variant_id             BIGINT PRIMARY KEY,
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
            PoN_AltDepth           INTEGER
        )
    """
    connection.execute(sql)

def drop_indexes(connection):
    log.logit("Dropping existing indexes on the variants table")
    sql = "DROP INDEX IF EXISTS chrom_pos_ref_alt_idx"
    connection.execute(sql)
    sql = "DROP INDEX IF EXISTS variant_id_idx"
    connection.execute(sql)
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

        log.logit("Generating the variant_id_idx")
        sql = """
            CREATE UNIQUE INDEX variant_id_idx
            ON variants ( variant_id )
        """
        connection.execute(sql)

        log.logit("Generating the batch_idx")
        sql = """
            CREATE INDEX batch_idx
            ON variants ( batch )
        """
        connection.execute(sql)

def setup_variants_table(connection):
    log.logit("Preparing the variant database file")
    ensure_variants_table(connection)
    drop_indexes(connection)
    return connection

def in_lsf_session():
    result = True if 'LSB_JOBID' in os.environ else False
    return result

def create_tmp_csv(work_dir, batch_number, chromosome):
    tmpdir = work_dir
    filename = None
    if in_lsf_session() and tmpdir == '/tmp':
        tmpdir = os.path.join('/tmp', f"{os.environ['LSB_JOBID']}.tmpdir")
    if chromosome is not None:
        filename = f"batch-{batch_number}.{chromosome}.csv.gz"
    else:
        filename = f"batch-{batch_number}.csv.gz"

    tmp_path = os.path.join(tmpdir, filename)
    if os.path.exists(tmp_path):
        log.logit(f"Deleting existing temporary csv file: {tmp_path}")
        os.unlink(tmp_path)

    return tmp_path

def get_variants(redis_db, batch_number, chromosome):
    set_name = f"batch:{batch_number}"
    scanner = None
    if chromosome:
        scanner = redis_db.sscan_iter(set_name, match=f"{chromosome}:*")
    else:
        scanner = redis_db.sscan_iter(set_name)

    for elem in scanner:
        variant_id = int(redis_db.get(elem))
        (chrom, pos, ref, alt) = elem.decode().split(':')
        pos = int(pos)
        start = pos
        stop = start + len(alt)
#        snp = 'TRUE' if len(ref) == len(alt) == 1 else 'FALSE'
        snp = True if len(ref) == len(alt) == 1 else False
        item = {
            'variant_id' : variant_id,
            'chrom'      : chrom,
            'pos'        : pos,
            'ref'        : ref,
            'alt'        : alt,
            'batch'      : batch_number,
            'start'      : start,
            'stop'       : stop,
            'snp'        : snp,
        }
        yield item

def csv_dump(redis_db, batch_number, chromosome, work_dir, debug):
    tmp_csv = create_tmp_csv(work_dir, batch_number, chromosome)
    log.logit(f"Placing variants into temporary csv file: {tmp_csv}")

    window_size = 100_000
    headers = ('variant_id', 'chrom', 'pos', 'ref', 'alt', 'batch', 'start', 'stop', 'snp')
    variants = get_variants(redis_db, batch_number, chromosome)
    count = 0
    with open(tmp_csv, 'wt') as fz, indent(4, quote=' >'):
        writer = csv.DictWriter(fz, fieldnames=headers)
        writer.writeheader()
        log.logit("Beginning to fetch variants into csv file")
        for variant in variants:
            count += 1
            if debug: log.logit(f"{count} -- {variant}")
            writer.writerow(variant)
            if count % window_size == 0:
                log.logit(f"# Variants Processed: {count}")

    log.logit(f"Finished fetching variants into csv file: {count} total variants")
    return (count, tmp_csv)

def duckdb_load_csv_file(duckdb_connection, temp_csv):
    sql = f"""
        COPY variants
        FROM read_csv(
            '{temp_csv}',
            delim=',',
            header=True,
            columns={{
                'variant_id' : 'BIGINT',
                'chrom'      : 'VARCHAR(5)',
                'pos'        : 'INTEGER',
                'ref'        : 'VARCHAR(255)',
                'alt'        : 'VARCHAR(255)',
                'batch'      : 'INTEGER',
                'start'      : 'INTEGER',
                'stop'       : 'INTEGER',
                'snp'        : 'BOOLEAN'
            }},
            compression=none
        )
    """
    log.logit(f"Starting to load csv into duckdb")
    duckdb_connection.execute(sql)
    log.logit(f"Finished loading csv into duckdb")

def insert_variants(redis_db, duckdb_connection, batch_number, chromosome, work_dir, debug):
    log.logit(f"Starting to ingest variants in batch: {batch_number}")
    (counts, temp_csv) = csv_dump(redis_db, batch_number, chromosome, work_dir, debug)
    duckdb_load_csv_file(duckdb_connection, temp_csv)
    return counts

def bulk_insert(duckdb_connection, entries):
    sql = "insert into variants ( variant_id, chrom, pos, ref, alt, batch, start, stop, snp ) values ( ?, ?, ?, ?, ?, ?, ?, ?, ?)"
    log.logit(f"Starting to bulk insert variant batch into duckdb ( {len(entries)} items)")
    duckdb.executemany(sql, entries, duckdb_connection)
    log.logit(f"Finished DuckDB insertion")

def insert_first_attempt(redis_db, duckdb_connection, batch_number, debug):
    redis_set_key = f"batch:{batch_number}"
    total_variants = int(redis_db.scard(redis_set_key))
    log.logit(f"Total # of variants to insert in batch {batch_number}: {total_variants}")
    log.logit("Inserting batch variants into duckdb variants table")
    window_size = 25_000
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
            if i % window_size == 0 and i != 0:
                bulk_insert(duckdb_connection, variant_batch)
                log.logit(f"# Variants Processed: {i}")
                variant_batch = []

        if len(variant_batch) > 0:
            bulk_insert(duckdb_connection, variant_batch)

    return { 'counts': total_variants }

def simple_bulk_insert(redis_db, duckdb_connection, batch_number, chromosome, window_size, debug):
    redis_set_key = f"batch:{batch_number}"
    log.logit(f"Inserting variants into duckdb variants table in increments of {window_size}")
    variants = get_variants(redis_db, batch_number, chromosome)
    attributes = ('variant_id', 'chrom', 'pos', 'ref', 'alt', 'batch', 'start', 'stop', 'snp')
    window = []
    count = 0
    with indent(4, quote=' >'):
        for (i, variant) in enumerate(variants):
            row = tuple(variant[attr] for attr in attributes)
            window.append( row )
            if debug: log.logit(f"{i} -- {row}")
            if i % window_size == 0 and i != 0:
                bulk_insert(duckdb_connection, window)
                log.logit(f"# Variants Processed: {i}")
                window = []
            count += 1

        if len(window) > 0:
            bulk_insert(duckdb_connection, window)

    return count

def ingest_variant_batch(duckdb_file, redis_host, redis_port, batch_number, chromosome, clobber, work_dir, window_size, debug):
    if chromosome:
        log.logit(f"Ingesting variants from redis batch: {batch_number} and chromosome: {chromosome} into {duckdb_file}", color="green")
    else:
        log.logit(f"Ingesting variants from redis batch: {batch_number} into {duckdb_file}", color="green")
    redis_db = db.redis_connect(redis_host, redis_port)
    duckdb_connection = db.duckdb_connect_rw(duckdb_file, clobber)
    setup_variants_table(duckdb_connection)
#    counts = insert_variants(redis_db, duckdb_connection, batch_number, chromosome, work_dir, debug)
    counts = simple_bulk_insert(redis_db, duckdb_connection, batch_number, chromosome, window_size, debug)
    create_indexes(duckdb_connection)
    duckdb_connection.close()
    log.logit(f"Finished ingesting variants")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def get_variants_from_table(duckdb_connection, batch_number):
    sql = f"SELECT chrom, pos, ref, alt FROM variants WHERE batch = {batch_number}"
    return duckdb_connection.sql(sql)

def dump_variant_batch(duckdb_file, header, batch_number, chromosome, work_dir, debug):
    import chip.vdbtools.importers.vcf as vcf
    log.logit(f"Dumping variants from batch: {batch_number} and chromosome: {chromosome} into a VCF file", color="green")
    duckdb_connection = db.duckdb_connect(duckdb_file)
    variants = get_variants_from_table(duckdb_connection, batch_number)
    vcf.write_variants_to_vcf(variants, header, batch_number, chromosome)
    duckdb_connection.close()
    log.logit(f"Finished dumping variants into VCF file")

#def merge_variant_tables():
    #  CREATE TABLE IF NOT EXISTS variants(
    #      variant_id             BIGINT PRIMARY KEY,
    #      chrom                  VARCHAR(5),
    #      pos                    INTEGER,
    #      ref                    VARCHAR(255) NOT NULL,
    #      alt                    VARCHAR(255) NOT NULL,
    #      snp                    BOOLEAN,
    #      qc_pass                BOOLEAN,
    #      batch                  INTEGER,
    #      start                  INTEGER,
    #      stop                   INTEGER
    #  );
    # ATTACH 'batch1.chr21.db' as chr21;
    # insert into variants select * from chr21.variants;
    # ATTACH 'batch1.chr22.db' as chr22;
    # insert into variants select * from chr22.variants;
    # create_indexes(connection)
