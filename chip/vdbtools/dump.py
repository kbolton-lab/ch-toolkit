def dump_variant_batch(duckdb_file, header, batch_number, chromosome, work_dir, debug):
    import chip.vdbtools.importers.variants as variants
    variants.dump_variant_batch(duckdb_file, header, batch_number, chromosome, work_dir, debug)