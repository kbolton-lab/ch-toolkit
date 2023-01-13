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
