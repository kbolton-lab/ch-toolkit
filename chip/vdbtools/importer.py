import sys

def import_vcf(db_path, input_vcf, caller, chromosome, variant_duckdb, sample_duckdb, clobber):
    dispatch = {
        'mutect'  : import_mutect,
        'vardict' : import_vardict
    }
    function = dispatch[caller]
    function(db_path, input_vcf, chromosome, variant_duckdb, sample_duckdb, clobber)

def import_mutect(db_path, input_vcf, chromosome, variant_duckdb, sample_duckdb, clobber):
    import chip.vdbtools.importers.callers as callers
    callers.insert_mutect_caller(db_path, input_vcf, chromosome, variant_duckdb, sample_duckdb, clobber)

def import_vardict(db_path, input_vcf, chromosome, variant_duckdb, sample_duckdb, clobber):
    import chip.vdbtools.importers.callers as callers
    callers.insert_vardict_caller(db_path, input_vcf, chromosome, variant_duckdb, sample_duckdb, clobber)

def import_variant_batch(duckdb_file, redis_host, redis_port, batch_number, chromosome, clobber, work_dir, window_size, debug):
    import chip.vdbtools.importers.variants as variants
    variants.ingest_variant_batch(duckdb_file, redis_host, redis_port, batch_number, chromosome, clobber, work_dir, window_size, debug)

def import_pon_pileup(variant_duckdb, pon_pileup, batch_number, chromosome, window_size, debug):
    import chip.vdbtools.importers.variants as variants
    variants.import_pon_pileup(variant_duckdb, pon_pileup, batch_number, chromosome, window_size, debug)
