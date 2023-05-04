import sys

def import_samples(samples_file, sample_duckdb, debug, clobber):
    import chip.vdbtools.importers.samples as samples
    samples.insert_samples(samples_file, sample_duckdb, debug, clobber)

def import_sample_variants(input_vcf, variant_db, batch_number, debug, clobber):
    import chip.vdbtools.importers.variants as variants
    variants.import_sample_variants(input_vcf, variant_db, batch_number, debug, clobber)

def import_vcf(db_path, input_vcf, caller, variant_db, sample_db, batch_number, clobber, debug):
    dispatch = {
        'mutect'  : import_mutect,
        'vardict' : import_vardict
    }
    function = dispatch[caller]
    function(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug)

def import_mutect(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug):
    import chip.vdbtools.importers.callers as callers
    callers.insert_mutect_caller(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug)

def import_vardict(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug):
    import chip.vdbtools.importers.callers as callers
    callers.insert_vardict_caller(db_path, input_vcf, variant_db, sample_db, batch_number, clobber, debug)

def import_caller_batch(db_path, caller_db, caller, batch_number, debug, clobber):
    import chip.vdbtools.importers.callers as callers
    callers.ingest_caller_batch(db_path, caller_db, caller, batch_number, debug, clobber)

def import_variant_batch(db_path, variant_db, batch_number, debug, clobber):
    import chip.vdbtools.importers.variants as variants
    variants.ingest_variant_batch(db_path, variant_db, batch_number, debug, clobber)

def import_pon_pileup(variant_duckdb, pon_pileup, batch_number, debug, clobber):
    import chip.vdbtools.importers.variants as variants
    variants.import_pon_pileup(variant_duckdb, pon_pileup, batch_number, debug, clobber)

def import_vep(annotation_db, variant_db, vep, batch_number, debug, clobber):
    import chip.vdbtools.importers.annotations as annotate
    annotate.import_vep(annotation_db, variant_db, vep, batch_number, debug, clobber)

def import_annotate_pd(annotation_db, annotate_pd, batch_number, debug):
    import chip.vdbtools.importers.annotations as annotate
    annotate.import_annotate_pd(annotation_db, annotate_pd, batch_number, debug)
