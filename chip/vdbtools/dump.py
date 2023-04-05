import sys

def dump_variant_batch(variant_db, header_type, batch_number, chromosome, debug):
    import chip.vdbtools.importers.variants as variants
    variants.dump_variant_batch(variant_db, header_type, batch_number, chromosome, debug)
