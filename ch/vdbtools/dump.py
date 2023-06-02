import sys

def dump_variant_batch(variant_db, header_type, batch_number, chromosome, debug):
    import ch.vdbtools.importers.variants as variants
    variants.dump_variant_batch(variant_db, header_type, batch_number, chromosome, debug)

def dump_variants_for_annotate_pd(annotation_db, batch_number, debug):
    import ch.vdbtools.importers.annotations as annotate
    annotate.dump_variants_batch(annotation_db, batch_number, debug)

def dump_ch_variants(mutect_db, vardict_db, annotation_db, debug):
    import ch.vdbtools.analysis.ch as ch
    ch.dump_ch_variants(mutect_db, vardict_db, annotation_db, debug)
