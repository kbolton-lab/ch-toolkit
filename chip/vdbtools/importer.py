import sys

def import_vcf(db_path, input_vcf, caller):
    dispatch = {
        'mutect'  : import_mutect,
        'vardict' : import_vardict
    }
    function = dispatch[caller]
    function(db_path, input_vcf)

def import_mutect(db_path, input_vcf):
    sys.exit("[err] Please implement me -- import_mutect !")

def import_vardict(db_path, input_vcf):
    sys.exit("[err] Please implement me -- import_vardict !")

def import_variant_batch(duckdb_file, redis_host, redis_port, batch_number, chromosome, clobber, work_dir, window_size, debug):
    import chip.vdbtools.importers.variants as variants
    variants.ingest_variant_batch(duckdb_file, redis_host, redis_port, batch_number, chromosome, clobber, work_dir, window_size, debug)
