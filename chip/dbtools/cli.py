import signal, importlib.resources
import click

from chip.version import __version__

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    '''A collection of db related tools for handling sample data.'''
    # to make this script/module behave nicely with unix pipes
    # http://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@cli.command('init-db', short_help="initialize a SNP/INDEL and SV candidate database")
@click.option('--db', 'db_path', default=None, type=click.Path(), required=True,
              help="Path to create sqlite database")
def init_db(db_path):
    with importlib.resources.path('chip.resources.sql', 'schema.sql') as f:
        schema_file = str(f.resolve())
        click.echo(schema_file)
