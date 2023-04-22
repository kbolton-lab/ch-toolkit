import os, sys, signal
import click
from clint.textui import puts, colored

from chip.version import __version__
import chip.utils.logger as log

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    '''A collection of db related tools for handling sample data.'''
    # to make this script/module behave nicely with unix pipes
    # http://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@cli.command('init', short_help="initialize a SNP/INDEL and SV candidate database")
@click.option('--db', 'db_path', default=None, type=click.Path(), required=True,
              help="Path to create sqlite database")
@click.option('--sample-name', 'sample_name', type=click.STRING, required=True,
              help="Sample Name for the variant database")
def init_db(db_path, sample_name):
    if os.path.exists(db_path):
        sys.exit(f"[err] Variant Database '{db_path}' already exists on filesystem!")
    import chip.dbtools.initialize as init
    init.create_db(db_path, sample_name)
    puts(colored.green(f"---> Successfully created variant database ({sample_name}): {db_path}"))

# Done
@cli.command('import-samples', short_help="Loads a CSV containing samples into samples database")
@click.option('--samples', '-s', type=click.Path(exists=True), required=True, help="A CSV file with the samples")
@click.option('--sdb', 'sample_duckdb', type=click.Path(), required=True, help="The duckdb database to store sample information")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def import_samples(samples, sample_duckdb, debug, clobber):
    """
    Loading samples into duckdb
    """
    import chip.vdbtools.importer as importer
    importer.import_samples(samples, sample_duckdb, debug, clobber)
    puts(colored.green(f"---> Successfully imported ({samples}) into {sample_duckdb}"))

# Done
@cli.command('import-sample-vcf', short_help="import a vcf file into sample variant database")
@click.option('--caller', 'caller',
              type=click.Choice(['mutect', 'vardict'], case_sensitive=False),
              required=True,
              help="Type of VCF file to import")
@click.option('--input-vcf', 'input_vcf', type=click.Path(exists=True), required=True, help="The VCF to be imported into the database")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
@click.option('--cdb', '-i', 'database', type=click.Path(), required=True, help="The duckdb database to import the caller data")
@click.option('--vdb', 'variant_db', type=click.Path(), required=True, help="The duckdb database to fetch variant ID from")
@click.option('--sdb', 'sample_db', type=click.Path(), required=True, help="The duckdb database to fetch sample ID from")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def import_vcf(caller, input_vcf, database, variant_db, sample_db, batch_number, clobber, debug):
    """
    variantdb is a path to a sample variant sqlite database.
    """
    import chip.vdbtools.importer as importer
    importer.import_vcf(database, input_vcf, caller, variant_db, sample_db, batch_number, clobber, debug)
    puts(colored.green(f"---> Successfully imported ({input_vcf}) into {database}"))

@cli.command('merge-batch-vcf', short_help="Combines all sample vcfs databases into a single database")
@click.option('--db-path', '-p', type=click.Path(exists=True), required=True, help="The path to where all the databases for this batch is stored")
@click.option('--cdb', 'caller_db', type=click.Path(), required=True, help="The variant database to import the batch into")
@click.option('--caller', 'caller',
              type=click.Choice(['mutect', 'vardict'], case_sensitive=False),
              required=True,
              help="Type of VCF file to import")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def merge_batch_vcf(db_path, caller_db, caller, batch_number, debug, clobber):
    """
    Ingest the variants in a batch into main variants database
    """
    import chip.vdbtools.importer as importer
    importer.import_caller_batch(db_path, caller_db, caller, batch_number, debug, clobber)
    log.logit(f"---> Successfully imported variant batch ({batch_number}) into {caller_db}", color="green")

@cli.command('register-variants', short_help="register the variants in a vcf file into redis")
@click.option('--input-vcf', '-i', 'input_vcf', type=click.Path(exists=True), required=True,
              help="The VCF to be imported into redis")
@click.option('--redis-host', '-h', type=click.STRING, required=True,
              help="The hostname of the redis server")
@click.option('--redis-port', '-p', type=click.IntRange(min=8000, max=8999), required=True,
              help="The port of the redis server")
@click.option('--batch-number', '-b', type=click.INT, required=True,
              help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True,
              help="Print extra debugging output")
def register_variants(input_vcf, redis_host, redis_port, batch_number, debug):
    """
    Registering variants into a redis database.
    """
    import chip.vdbtools.register as register
    register.import_vcf(input_vcf, redis_host, redis_port, batch_number, debug)
    log.logit(f"---> Successfully imported ({input_vcf}) into redis://{redis_host}:{redis_port}", color="green")

# Done - Same as: register-variants
@cli.command('register-sample-variants', short_help="Register the variants for a VCF file into a variant database")
@click.option('--input-vcf', '-i', 'input_vcf', type=click.Path(exists=True), required=True, help="The VCF to be imported")
@click.option('--db', 'duck_db', type=click.Path(), required=True, help="The variant database to import the VCF into")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def register_sample_variants(input_vcf, duck_db, batch_number, debug, clobber):
    """
    Registering variants into the duckdb database.
    """
    import chip.vdbtools.importer as importer
    importer.import_sample_variants(input_vcf, duck_db, batch_number, debug, clobber)
    log.logit(f"---> Successfully imported ({input_vcf})", color="green")

# Done - Same as: ingest-variants
@cli.command('merge-batch-variants', short_help="Combines all sample variant databases into a single database")
@click.option('--db-path', '-p', type=click.Path(exists=True), required=True, help="The path to where all the databases for this batch is stored")
@click.option('--vdb', 'variant_db', type=click.Path(), required=True, help="The variant database to import the batch into")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def merge_batch_variants(db_path, variant_db, batch_number, debug, clobber):
    """
    Ingest the variants in a batch into main variants database
    """
    import chip.vdbtools.importer as importer
    importer.import_variant_batch_(db_path, variant_db, batch_number, debug, clobber)
    log.logit(f"---> Successfully imported variant batch ({batch_number}) into {variant_db}", color="green")

@cli.command('ingest-variants', short_help="ingest the variants in a batch from redis into duckdb")
@click.option('--db', '-i', 'database', type=click.Path(), required=True,
              help="The duckdb database to import the variant set from redis")
@click.option('--redis-host', '-h', type=click.STRING, required=True,
              help="The hostname of the redis server")
@click.option('--redis-port', '-p', type=click.IntRange(min=8000, max=8999), required=True,
              help="The port of the redis server")
@click.option('--batch-number', '-b', type=click.INT, required=True,
              help="The batch number of this import set")
@click.option('--chromosome', '-c', type=click.STRING, default=None,
              help="The chromosome set of interest")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True,
              help="If exists, delete existing duckdb file and then start from scratch")
@click.option('--work-dir', '-d', type=click.Path(exists=True), default="/tmp", show_default=True, required=False,
              help="The the working directory to create temporary files (usually the OS temp directory)")
@click.option('--window-size', '-w', type=click.INT, default=10_000, show_default=True, required=False,
              help="The variant window size when bulk ingesting variants by executemany")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True,
              help="Print extra debugging output")
def ingest_variants(database, redis_host, redis_port, batch_number, chromosome, clobber, work_dir, window_size, debug):
    """
    Ingest the variants in a batch from redis into duckdb's variant table
    """
    import chip.vdbtools.importer as importer
    importer.import_variant_batch(database, redis_host, redis_port, batch_number, chromosome, clobber, work_dir, window_size, debug)
    log.logit(f"---> Successfully imported variant batch ({batch_number}) into {database}", color="green")

# Done
@cli.command('dump-variants', short_help="dumps all variants inside duckdb into a VCF file")
@click.option('--vdb', 'variant_db', type=click.Path(exists=True), required=True, help="The duckdb database to dump the variants from")
@click.option('--header-type', '-t', type=click.Choice(['simple', 'dummy'], case_sensitive=False), required=True, default="dummy",
                                    help="A pre-existing header type e.g. simple, mutect, vardict, complex, etc.")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--chromosome', '-c', type=click.STRING, default=None, help="The chromosome set of interest")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def dump_variants(variant_db, header_type, batch_number, chromosome, debug):
    """
    Dumps the variants from duckdb into a VCF file
    """
    import chip.vdbtools.dump as dump
    dump.dump_variant_batch(variant_db, header_type, batch_number, chromosome, debug)
    log.logit(f"---> Successfully dumped variant batch ({batch_number}) from {variant_db}", color="green")

# Done
@cli.command('import-pon-pileup', short_help="updates variants inside duckdb with PoN pileup information")
@click.option('--vdb', 'variant_db', type=click.Path(exists=True), required=True, help="The duckdb database to fetch variant ID from")
@click.option('--pon-pileup', '-p', 'pon_pileup', type=click.Path(exists=True), required=True, help="The pon pileup VCF to be imported into the variant database")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def import_pon_pileup(variant_db, pon_pileup, batch_number, debug, clobber):
    """
    Dumps the panel of normal pileup information from into variants duckdb
    """
    import chip.vdbtools.importer as importer
    importer.import_pon_pileup(variant_db, pon_pileup, batch_number, debug, clobber)
    log.logit(f"---> Successfully imported PoN Pileup from batch ({batch_number}) into {variant_db}", color="green")

# Done
@cli.command('calculate-fishers-test', short_help="Updates the variants inside Mutect or Vardict tables with p-value from Fisher's Exact Test")
@click.option('--vdb', 'variant_db', type=click.Path(exists=True), required=True, help="The duckdb database to fetch variant PoN Ref Depth and Alt Depth from")
@click.option('--cdb', 'caller_db', type=click.Path(exists=True), required=True, help="The duckdb database to fetch variant caller information from")
@click.option('--caller', 'caller',
              type=click.Choice(['mutect', 'vardict'], case_sensitive=False),
              required=True,
              help="Type of VCF file to import")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def calculate_fishers_test(variant_db, caller_db, caller, batch_number, debug):
    """
    Calculates the Fisher's Exact Test for all Variants within the Variant Caller duckdb
    """
    import chip.vdbtools.importers.callers as callers
    callers.annotate_fisher_test(variant_db, caller_db, caller, batch_number, debug)
    log.logit(f"---> Successfully calculated the Fisher's Exact Test for variants within ({batch_number}) and {caller_db}", color="green")

# TODO
@cli.command('import-vep', short_help="updates variants inside duckdb with VEP information")
@click.option('--adb', 'annotation_db', type=click.Path(), required=True, help="The duckdb database to store the annotation information")
@click.option('--vdb', 'variant_db', type=click.Path(exists=True), required=True, help="The duckdb database to fetch variant key from")
@click.option('--vep', '-v', type=click.Path(exists=True), required=True, help="The VEP TSV to be imported into the annotation database")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def import_vep(annotation_db, variant_db, vep, batch_number, debug, clobber):
    """
    Dumps the vep information into an annotation duckdb
    """
    import chip.vdbtools.importer as importer
    importer.import_vep(annotation_db, variant_db, vep, batch_number, debug, clobber)
    log.logit(f"---> Successfully imported VEP from batch ({batch_number}) into {annotation_db}", color="green")

@cli.command('dump-annotations', short_help="dumps all variant annotations inside duckdb into a CSV file")
@click.option('--adb', 'annotation_db', type=click.Path(exists=True), required=True, help="The duckdb database to dump the variant annotations from")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def dump_variants_for_annotate_pd(annotation_db, batch_number, debug):
    """
    Dumps the variant annotations from duckdb into a CSV file
    """
    import chip.vdbtools.dump as dump
    dump.dump_variants_for_annotate_pd(annotation_db, batch_number, debug)
    log.logit(f"---> Successfully dumped variant annotations batch ({batch_number}) from {annotation_db}", color="green")

# TODO
@cli.command('import-annotate-pd', short_help="annotates variants with their pathogenicity")
@click.option('--adb', 'annotation_db', type=click.Path(exists=True), required=True, help="The duckdb database to store the annotation information")
@click.option('--pd', '-p', 'annotate_pd', type=click.Path(exists=True), required=True, help="The CSV File produced from AnnotatePD to be imported into the annotation database")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def import_annotate_pd(annotation_db, annotate_pd, batch_number, debug):
    """
    Dumps the vep information into an annotation duckdb
    """
    import chip.vdbtools.importer as importer
    importer.import_annotate_pd(annotation_db, annotate_pd, batch_number, debug)
    log.logit(f"---> Successfully annotated variants from batch ({batch_number}) in {annotation_db}", color="green")
