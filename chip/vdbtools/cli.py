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

@cli.command('import-vcf', short_help="import a vcf file into sample variant database")
@click.option('--caller', 'caller',
              type=click.Choice(['lofreq', 'mutect', 'vardict', 'pindel'], case_sensitive=False),
              required=True,
              help="Type of VCF file to import")
@click.option('--input-vcf', 'input_vcf', type=click.Path(exists=True), required=True,
              help="The VCF to be imported into the database")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True,
              help="If exists, delete existing duckdb file and then start from scratch")
@click.option('--db', '-i', 'database', type=click.Path(), required=True,
              help="The duckdb database to import the caller data")
@click.option('--vdb', 'variant_db', type=click.Path(), required=True,
              help="The duckdb database to fetch variant ID from")
@click.option('--sdb', 'sample_db', type=click.Path(), required=True,
              help="The duckdb database to fetch sample ID from")
@click.option('--window-size', '-w', type=click.INT, default=10_000, show_default=True, required=False,
              help="The variant window size when bulk updating variants by executemany")
@click.option('--chromosome', '-c', type=click.STRING, default=None,
              help="The chromosome set of interest")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True,
              help="Print extra debugging output")
def import_vcf(caller, input_vcf, database, chromosome, variant_duckdb, sample_duckdb, clobber, debug, window_size):
    """
    variantdb is a path to a sample variant sqlite database.
    """
    import chip.vdbtools.importer as importer
    importer.import_vcf(database, input_vcf, caller, chromosome, variant_duckdb, sample_duckdb, clobber, debug, window_size)
    puts(colored.green(f"---> Successfully imported ({input_vcf}) into {database}"))

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

## TESTING
@cli.command('register-sample-variants', short_help="Register the variants for a VCF file into a variant database")
@click.option('--input-vcf', '-i', 'input_vcf', type=click.Path(exists=True), required=True, help="The VCF to be imported")
@click.option('--db', 'duck_db', type=click.Path(), required=True, help="The variant database to import the VCF into")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
@click.option('--work-dir', '-w', type=click.Path(exists=True), default="/tmp", show_default=True, required=False,
              help="The the working directory to create temporary files (usually the OS temp directory)")
def register_sample_variants(input_vcf, duck_db, batch_number, debug, clobber, work_dir):
    """
    Registering variants into the duckdb database.
    """
    import chip.vdbtools.importer as importer
    importer.import_sample_variants(input_vcf, duck_db, batch_number, debug, clobber, work_dir)
    log.logit(f"---> Successfully imported ({input_vcf})", color="green")

@cli.command('ingest-batch-variants', short_help="Combines all sample variant databases into a single database")
@click.option('--db-path', '-p', type=click.Path(exists=True), required=True, help="The path to where all the databases for this batch is stored")
@click.option('--vdb', 'variant_db', type=click.Path(), required=True, help="The variant database to import the batch into")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
@click.option('--work-dir', '-w', type=click.Path(exists=True), default="/tmp", show_default=True, required=False,
              help="The the working directory to create temporary files (usually the OS temp directory)")
def ingest_batch_variants(db_path, variant_db, batch_number, debug, clobber, work_dir):
    """
    Ingest the variants in a batch into main variants database
    """
    import chip.vdbtools.importer as importer
    importer.import_variant_batch_(db_path, variant_db, batch_number, debug, clobber, work_dir)
    log.logit(f"---> Successfully imported variant batch ({batch_number}) into {variant_db}", color="green")
###

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

@cli.command('dump-variants', short_help="dumps all variants inside duckdb into a VCF file")
@click.option('--vdb', 'variant_db', type=click.Path(exists=True), required=True, help="The duckdb database to dump the variants from")
@click.option('--header-type', '-t', type=click.STRING, default=None, required=True, help="A pre-existing header type e.g. simple, mutect, vardict, complex, etc.")
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

@cli.command('import-pon-pileup', short_help="updates variants inside duckdb with PoN pileup information")
@click.option('--vdb', 'variant_duckdb', type=click.Path(exists=True), required=True,
              help="The duckdb database to fetch variant ID from")
@click.option('--pon-pileup', '-p', 'pon_pileup', type=click.Path(exists=True), required=True,
              help="The pon pileup VCF to be imported into the variant database")
@click.option('--batch-number', '-b', type=click.INT, required=True,
              help="The batch number of this variant set")
@click.option('--chromosome', '-c', type=click.STRING, default=None,
              help="The chromosome set of interest")
@click.option('--window-size', '-w', type=click.INT, default=10_000, show_default=True, required=False,
              help="The variant window size when bulk updating variants by executemany")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True,
              help="Print extra debugging output")
def import_pon_pileup(variant_duckdb, pon_pileup, batch_number, chromosome, window_size, debug):
    """
    Dumps the panel of normal pileup information from into variants duckdb
    """
    import chip.vdbtools.importer as importer
    importer.import_pon_pileup(variant_duckdb, pon_pileup, batch_number, chromosome, window_size, debug)
    log.logit(f"---> Successfully imported PoN Pileup from batch ({batch_number}) into {variant_duckdb}", color="green")

@cli.command('calculate-fishers-test', short_help="Updates the variants inside Mutect or Vardict tables with p-value from Fisher's Exact Test")
@click.option('--vdb', 'variant_duckdb', type=click.Path(exists=True), required=True,
              help="The duckdb database to fetch variant PoN Ref Depth and Alt Depth from")
@click.option('--cdb', 'caller_duckdb', type=click.Path(exists=True), required=True,
              help="The duckdb database to fetch variant caller information from")
@click.option('--caller', 'caller',
              type=click.Choice(['lofreq', 'mutect', 'vardict', 'pindel'], case_sensitive=False),
              required=True,
              help="Type of VCF file to import")
@click.option('--batch-number', '-b', type=click.INT, required=True,
              help="The batch number of this variant set")
@click.option('--chromosome', '-c', type=click.STRING, default=None,
              help="The chromosome set of interest")
@click.option('--window-size', '-w', type=click.INT, default=10_000, show_default=True, required=False,
              help="The variant window size when bulk updating variants by executemany")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True,
              help="Print extra debugging output")
def calculate_fishers_test(variant_duckdb, caller_duckdb, caller, batch_number, chromosome, window_size, debug):
    """
    Calculates the Fisher's Exact Test for all Variants within the Variant Caller duckdb
    """
    import chip.vdbtools.importers.callers as callers
    callers.annotate_fisher_test(variant_duckdb, caller_duckdb, caller, batch_number, chromosome, window_size, debug)
    log.logit(f"---> Successfully calculated the Fisher's Exact Test for variants within ({batch_number}) and {caller_duckdb}", color="green")
