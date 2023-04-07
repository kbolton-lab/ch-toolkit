import sys

def import_samples(samples_file, sample_duckdb, debug, clobber):
    import chip.vdbtools.importers.samples as samples
    samples.insert_samples(samples_file, sample_duckdb, debug, clobber)

def import_sample_variants(input_vcf, duck_db, batch_number, debug, clobber):
    import chip.vdbtools.importers.variants as variants
    variants.import_sample_variants(input_vcf, duck_db, batch_number, debug, clobber)

def import_vcf(db_path, input_vcf, caller, variant_db, sample_db, batch_number, clobber, debug):
    dispatch = {
        'mutect'  : import_mutect,
        'vardict' : import_vardict
    }
    function = dispatch[caller]
    function(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug)

def import_mutect(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug):
    import chip.vdbtools.importers.callers as callers
    callers.insert_mutect_caller_(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug)

def import_vardict(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug):
    import chip.vdbtools.importers.callers as callers
    callers.insert_vardict_caller(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug)

def import_variant_batch(duckdb_file, redis_host, redis_port, batch_number, chromosome, clobber, work_dir, window_size, debug):
    import chip.vdbtools.importers.variants as variants
    variants.ingest_variant_batch(duckdb_file, redis_host, redis_port, batch_number, chromosome, clobber, work_dir, window_size, debug)

def import_variant_batch_(db_path, variant_db, batch_number, debug, clobber):
    import chip.vdbtools.importers.variants as variants
    variants.ingest_variant_batch_(db_path, variant_db, batch_number, debug, clobber)

def import_pon_pileup(variant_duckdb, pon_pileup, batch_number, debug, clobber):
    import chip.vdbtools.importers.variants as variants
    variants.import_pon_pileup(variant_duckdb, pon_pileup, batch_number, debug, clobber)
