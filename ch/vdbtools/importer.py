import sys

def import_samples(samples_file, sample_duckdb, batch_number, debug, clobber):
    import ch.vdbtools.handlers.samples as samples
    samples.insert_samples(samples_file, sample_duckdb, batch_number, debug, clobber)

def import_sample_variants(input_vcf, variant_db, batch_number, debug, clobber):
    import ch.vdbtools.handlers.variants as variants
    variants.import_sample_variants(input_vcf, variant_db, batch_number, debug, clobber)

def import_vcf(db_path, input_vcf, caller, batch_number, clobber, debug):
    dispatch = {
        'mutect'  : import_mutect,
        'vardict' : import_vardict
    }
    function = dispatch[caller]
    function(db_path, input_vcf, batch_number, clobber, debug)

def import_mutect(db_path, input_vcf, batch_number, clobber, debug):
    import ch.vdbtools.handlers.callers as callers
    callers.insert_mutect_caller(db_path, input_vcf, batch_number, clobber, debug)

def import_vardict(db_path, input_vcf, batch_number, clobber, debug):
    import ch.vdbtools.handlers.callers as callers
    callers.insert_vardict_caller(db_path, input_vcf, batch_number, clobber, debug)

def import_caller_batch(db_path, caller_db, variant_db, sample_db, caller, batch_number, debug, clobber):
    import ch.vdbtools.handlers.callers as callers
    callers.insert_caller_batch(db_path, caller_db, variant_db, sample_db, caller, batch_number, debug, clobber)

def import_variant_batch(db_path, variant_db, batch_number, debug, clobber):
    import ch.vdbtools.handlers.variants as variants
    variants.insert_variant_batch(db_path, variant_db, batch_number, debug, clobber)

def import_pon_pileup(pileup_db, variant_db, pon_pileup, batch_number, debug, clobber):
    import ch.vdbtools.handlers.variants as variants
    variants.import_pon_pileup(pileup_db, variant_db, pon_pileup, batch_number, debug, clobber)

def import_vep(annotation_db, variant_db, vep, batch_number, debug, clobber):
    import ch.vdbtools.handlers.annotations as annotate
    annotate.import_vep(annotation_db, variant_db, vep, batch_number, debug, clobber)

def import_annotate_pd(annotation_db, annotate_pd, batch_number, debug):
    import ch.vdbtools.handlers.annotations as annotate
    annotate.import_annotate_pd(annotation_db, annotate_pd, batch_number, debug)
