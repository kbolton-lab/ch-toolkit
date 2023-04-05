import os, csv, glob
import vcfpy
import duckdb

import chip.utils.logger as log
import chip.utils.database as db
import chip.utils.csv_utils as csv_utils
from clint.textui import indent, puts_err, puts

#            variant_id             BIGINT PRIMARY KEY,
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
            PoN_AltDepth           INTEGER
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

def setup_variants_table(connection):
    log.logit("Preparing the variant database file")
    ensure_variants_table(connection)
    drop_indexes(connection)
    return connection

def get_variants(redis_db, batch_number, chromosome):
    set_name = f"batch:{batch_number}"
    # scanner = None
    elems = [elem.decode() for elem in redis_db.smembers(set_name)]
    # if chromosome:
    #     scanner = redis_db.keys(f"{chromosome}:*")
    #     # scanner = redis_db.sscan_iter(set_name, match=f"{chromosome}:*")
    # else:
    #     scanner = redis_db.keys(f"chr*")
    if chromosome is not None:
        elems = [elem for elem in elems if chromosome in elem]
    pipe = redis_db.pipeline()
    for elem in elems:
        pipe.get(elem)
    variant_ids = pipe.execute()
    res = []
    for i, elem in enumerate(elems):
        variant_id = int(variant_ids[i])
        (chrom, pos, ref, alt) = elem.split(':')
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
            'snp'        : snp,
            'qc_pass'    : None,
            'batch'      : batch_number,
            'start'      : start,
            'stop'       : stop,
            'PoN_RefDepth': None,
            'PoN_AltDepth': None
        }
        res.append(item)
        #yield item
    return res

def csv_dump(redis_db, batch_number, chromosome, work_dir, debug):
    tmp_csv = csv_utils.create_tmp_csv(work_dir, batch_number, chromosome)
    log.logit(f"Placing variants into temporary csv file: {tmp_csv}")
    window_size = 100_000
    headers = ('variant_id', 'chrom', 'pos', 'ref', 'alt', 'snp', 'qc_pass', 'batch', 'start', 'stop', 'PoN_RefDepth', 'PoN_AltDepth')
    variants = get_variants(redis_db, batch_number, chromosome)
    count = len(variants)
    with open(tmp_csv, 'wt') as fz, indent(4, quote=' >'):
        writer = csv.DictWriter(fz, fieldnames=headers)
        writer.writeheader()
        log.logit("Beginning to fetch variants into csv file")
        writer.writerows(variants)
    log.logit(f"Finished fetching variants into csv file: {count} total variants")
    return (count, tmp_csv)

def insert_variants(redis_db, duckdb_connection, batch_number, chromosome, work_dir, debug):
    log.logit(f"Starting to ingest variants in batch: {batch_number}")
    (counts, temp_csv) = csv_dump(redis_db, batch_number, chromosome, work_dir, debug)
    csv_utils.duckdb_load_csv_file(duckdb_connection, temp_csv, "variants")
    return counts

# def bulk_insert(duckdb_connection, entries):
#     sql = "insert into variants ( variant_id, chrom, pos, ref, alt, batch, start, stop, snp ) values ( ?, ?, ?, ?, ?, ?, ?, ?, ?)"
#     log.logit(f"Starting to bulk insert variant batch into duckdb ( {len(entries)} items)")
#     duckdb.executemany(sql, entries, duckdb_connection)
#     log.logit(f"Finished DuckDB insertion")

def bulk_update_pileup(duckdb_connection, entries):
    sql = f"UPDATE variants SET PoN_RefDepth = ?, PoN_AltDepth = ? WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?"
    log.logit(f"Starting to bulk update variant batch into duckdb ( {len(entries)} items)")
    duckdb.executemany(sql, entries, duckdb_connection)
    log.logit(f"Finished Updating the Variants in DuckDB")

# def insert_first_attempt(redis_db, duckdb_connection, batch_number, debug):
#     redis_set_key = f"batch:{batch_number}"
#     total_variants = int(redis_db.scard(redis_set_key))
#     log.logit(f"Total # of variants to insert in batch {batch_number}: {total_variants}")
#     log.logit("Inserting batch variants into duckdb variants table")
#     window_size = 25_000
#     variant_batch = []
#     with indent(4, quote=' >'):
#         for (i, key) in enumerate(redis_db.sscan_iter(redis_set_key)):
#             variant_id = int(redis_db.get(key))
#             (chrom, pos, ref, alt) = key.decode().split(':')
#             pos = int(pos)
#             start = pos
#             stop = start + len(alt)
#             item = (variant_id, chrom, pos, ref, alt, batch_number, start, stop)
#             variant_batch.append( item )
#             if debug: log.logit(f"{i} -- {item}")
#             if i % window_size == 0 and i != 0:
#                 bulk_insert(duckdb_connection, variant_batch)
#                 log.logit(f"# Variants Processed: {i}")
#                 variant_batch = []
#
#         if len(variant_batch) > 0:
#             bulk_insert(duckdb_connection, variant_batch)
#     return { 'counts': total_variants }

# def simple_bulk_insert(redis_db, duckdb_connection, batch_number, chromosome, window_size, debug):
#     redis_set_key = f"batch:{batch_number}"
#     log.logit(f"Inserting variants into duckdb variants table in increments of {window_size}")
#     variants = get_variants(redis_db, batch_number, chromosome)
#     attributes = ('variant_id', 'chrom', 'pos', 'ref', 'alt', 'batch', 'start', 'stop', 'snp')
#     window = []
#     count = 0
#     with indent(4, quote=' >'):
#         for (i, variant) in enumerate(variants):
#             row = tuple(variant[attr] for attr in attributes)
#             window.append( row )
#             if debug: log.logit(f"{i} -- {row}")
#             if i % window_size == 0 and i != 0:
#                 bulk_insert(duckdb_connection, window)
#                 log.logit(f"# Variants Processed: {i}")
#                 window = []
#             count += 1
#
#         if len(window) > 0:
#             bulk_insert(duckdb_connection, window)
#     return count

def merge_variants_tables(db_path, connection, batch_number, debug, clobber):
    log.logit(f'Merging Sample Variants')
    with indent(4, quote=' >'):
        for i, file in enumerate(glob.glob(db_path + "*.db")):
            sample_name = os.path.basename(file).split('.')[1]
            log.logit(f"Merging: {file}")
            connection.execute(f"ATTACH \'{file}\' as sample_{i}")
            sql = f"""
                INSERT INTO variants SELECT s.*
                FROM sample_{i}.variants s
                WHERE NOT EXISTS (
                    SELECT *
                    FROM variants v
                    WHERE (
                        s.chrom = v.chrom
                        AND s.pos = v.pos
                        AND s.ref = v.ref
                        AND s.alt = v.alt
                    )
                )
            """
            connection.sql(sql)
    log.logit(f"Finished merging all tables from: {db_path}")
    #log.logit(f"Variants Processed - Total: {counts}", color="green")

    #     connection.execute(f"ATTACH 'batch{batch_number}.chr{i}.db' as chr{i}")
    #     connection.execute(f"insert into variants select * from chr{i}.variants")

def merge_variants_tables_chromosomes(variants_db, batch_number, clobber, debug):
    connection = db.duckdb_connect_rw(variants_db, clobber)
    setup_variants_table(connection)
    chromosomes = [str(c) for c in range(1,23)]
    chromosomes.extend(['X', 'Y'])
    for i in chromosomes:
        log.logit(f'Merging chromosome: {i}')
        connection.execute(f"ATTACH 'batch{batch_number}.chr{i}.db' as chr{i}")
        connection.execute(f"insert into variants select * from chr{i}.variants")
    log.logit("Finished merging chromosomes")
    log.logit("Creating Indexes")
    create_indexes(connection)
    connection.close()
    log.logit(f"Finished combining all variants into: {variants_db}")
    #log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def ingest_variant_batch(duckdb_file, redis_host, redis_port, batch_number, chromosome, clobber, work_dir, window_size, debug):
    if chromosome:
        log.logit(f"Ingesting variants from redis batch: {batch_number} and chromosome: {chromosome} into {duckdb_file}", color="green")
    else:
        log.logit(f"Ingesting variants from redis batch: {batch_number} into {duckdb_file}", color="green")
    redis_db = db.redis_connect(redis_host, redis_port)
    duckdb_connection = db.duckdb_connect_rw(duckdb_file, clobber)
    setup_variants_table(duckdb_connection)
    counts = insert_variants(redis_db, duckdb_connection, batch_number, chromosome, work_dir, debug)
    #counts = simple_bulk_insert(redis_db, duckdb_connection, batch_number, chromosome, window_size, debug)
    create_indexes(duckdb_connection)
    duckdb_connection.close()
    log.logit(f"Finished ingesting variants")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def ingest_variant_batch_(db_path, variant_db, batch_number, debug, clobber, work_dir):
    log.logit(f"Ingesting variants from batch: {batch_number} into {variant_db}", color="green")
    connection = db.duckdb_connect_rw(variant_db, clobber)
    setup_variants_table(connection)
    merge_variants_tables(db_path, connection, batch_number, debug, clobber)
    create_indexes(connection)
    connection.close()
    log.logit(f"Finished ingesting variants")
    #log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def import_sample_variants(input_vcf, duck_db, batch_number, debug, clobber, work_dir):
    log.logit(f"Registering variants from file: {input_vcf}", color="green")
    connection = db.duckdb_connect_rw(duck_db, clobber)
    setup_variants_table(connection)
    headers = ('chrom', 'pos', 'ref', 'alt', 'snp', 'qc_pass', 'batch', 'start', 'stop', 'PoN_RefDepth', 'PoN_AltDepth')
    counts, tmp_csv = csv_utils.vcf_to_csv(input_vcf, "variants", headers, batch_number, work_dir, debug)
    csv_utils.duckdb_load_csv_file(connection, tmp_csv, "variants")
    create_indexes(connection)
    connection.close()
    log.logit(f"Finished registering variants")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")

def get_variants_from_table(connection, batch_number, chromosome):
    if chromosome != None:
        log.logit(f"Grabbing variants from batch: {batch_number} and chromosome: {chromosome} from the database")
        sql = f"SELECT variant_id, chrom, pos, ref, alt FROM variants WHERE batch = {batch_number} AND chrom = \'{chromosome}\'"
    else:
        log.logit(f"Grabbing variants from batch: {batch_number} from the database")
        sql = f"SELECT ROWID, chrom, pos, ref, alt FROM variants WHERE batch = {batch_number}"
    return connection.sql(sql)

def dump_variant_batch(variant_db, header, batch_number, chromosome, debug):
    import chip.vdbtools.importers.vcf as vcf
    if chromosome is None:
        log.logit(f"Dumping batch: {batch_number} variants from: {variant_db} into a VCF file", color="green")
    else:
        log.logit(f"Dumping batch: {batch_number} and chromosome: {chromosome} variants from: {variant_db} into a VCF file", color="green")
    connection = db.duckdb_connect(variant_db)
    variants = get_variants_from_table(connection, batch_number, chromosome)
    vcf.write_variants_to_vcf(variants, header, batch_number, chromosome, debug)
    connection.close()
    log.logit(f"Finished dumping variants into VCF file")
    log.logit(f"All Done!", color="green")

def annotate_variants_with_pon_pileup(connection, pon_pileup, chromosome, window_size, debug):
    log.logit(f"Importing pileup information in increments of {window_size}")
    reader = vcfpy.Reader.from_path(pon_pileup)
    variants = reader.fetch(chromosome) if chromosome != None else reader
    attributes = ('PoN_RefDepth', 'PoN_AltDepth')
    window = []
    count = 0
    with indent(4, quote=' >'):
        for (i, variant) in enumerate(variants):
            row = (variant.INFO.get('PON_RefDepth')[0], variant.INFO.get('PON_AltDepth')[0], variant.CHROM, variant.POS, variant.REF, variant.ALT[0].value)
            if debug: log.logit(f"{i} -- {row}")
            window.append( row )
            if i % window_size == 0 and i != 0:
                bulk_update_pileup(connection, window)
                log.logit(f"# Variants Processed: {i}")
                window = []
            count += 1
        if len(window) > 0:
            bulk_update_pileup(connection, window)
    return count

def import_pon_pileup(variant_duckdb, pon_pileup, batch_number, chromosome, window_size, debug):
    if chromosome:
        log.logit(f"Adding pileup from batch: {batch_number} and chromosome: {chromosome} into {variant_duckdb}", color="green")
    else:
        log.logit(f"Adding pileup from batch: {batch_number} into {variant_duckdb}", color="green")
    variant_duckdb_connection = db.duckdb_connect_rw(variant_duckdb, False)
    counts = annotate_variants_with_pon_pileup(variant_duckdb_connection, pon_pileup, chromosome, window_size, debug)
    variant_duckdb_connection.close()
    log.logit(f"Finished importing PoN pileup information")
    log.logit(f"Variants Processed - Total: {counts}", color="green")
    log.logit(f"All Done!", color="green")
