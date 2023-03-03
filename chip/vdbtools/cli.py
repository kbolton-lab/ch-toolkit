import os, sys, signal
import click
from clint.textui import puts, colored

from chip.version import __version__

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

@cli.command('import-vcf', short_help="import a vcf file into sample variant database")
@click.option('--caller', 'caller',
              type=click.Choice(['lofreq', 'mutect', 'vardict', 'pindel'], case_sensitive=False),
              required=True,
              help="Type of VCF file to import")
@click.option('--input-vcf', 'input_vcf', type=click.Path(exists=True), required=True,
              help="The VCF to be imported into the database")
@click.argument('variantdb', type=click.Path(exists=True), nargs=1)
def import_vcf(caller, input_vcf, variantdb):
    """
    variantdb is a path to a sample variant sqlite database.
    """
    import chip.vdbtools.importer as importer
    importer.import_vcf(variantdb, input_vcf, caller)
    puts(colored.green(f"---> Successfully imported ({input_vcf}) into {variantdb}"))

@cli.command('register-variants', short_help="register the variants in vcf file into redis")
@click.option('--input-vcf', 'input_vcf', type=click.Path(exists=True), required=True,
              help="The VCF to be imported into redis")
@click.option('--redis-host', type=click.STRING, required=True,
              help="The hostname of the redis server")
@click.option('--redis-port', type=click.IntRange(min=8000, max=8999), required=True,
              help="The port of the redis server")
@click.option('--batch-number', type=click.INT, required=True,
              help="The batch number of this import set")
def register_variants(input_vcf, redis_host, redis_port, batch_number):
    """
    variantdb is a path to a sample variant sqlite database.
    """
    import chip.vdbtools.register as register
    register.import_vcf(input_vcf, redis_host, redis_port, batch_number)
    puts(colored.green(f"---> Successfully imported ({input_vcf}) into {redis_host}"))
