import os, sys, signal
import click
from clint.textui import puts, colored

from ch.version import __version__
import ch.utils.logger as log

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    '''A collection of db related tools for handling sample data.'''
    # to make this script/module behave nicely with unix pipes
    # http://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@cli.command('import-samples', short_help="Loads a CSV containing samples into samples database")
@click.option('--samples', '-s', type=click.Path(exists=True), required=True, help="A CSV file with the samples")
@click.option('--sdb', 'sample_duckdb', type=click.Path(), required=True, help="The duckdb database to store sample information")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def import_samples(samples, sample_duckdb, batch_number, debug, clobber):
    """
    Loading samples into duckdb
    """
    import ch.vdbtools.importer as importer
    importer.import_samples(samples, sample_duckdb, batch_number, debug, clobber)
    puts(colored.green(f"---> Successfully imported ({samples}) from batch {batch_number} into {sample_duckdb}"))

@cli.command('import-sample-vcf', short_help="import a vcf file into sample variant database")
@click.option('--caller', 'caller',
              type=click.Choice(['mutect', 'vardict'], case_sensitive=False),
              required=True,
              help="Type of VCF file to import")
@click.option('--input-vcf', 'input_vcf', type=click.Path(exists=True), required=True, help="The VCF to be imported into the database")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
@click.option('--cdb', '-i', 'database', type=click.Path(), required=True, help="The duckdb database to import the caller data")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def import_vcf(caller, input_vcf, database, batch_number, clobber, debug):
    """
    variantdb is a path to a sample variant sqlite database.
    """
    import ch.vdbtools.importer as importer
    importer.import_vcf(database, input_vcf, caller, batch_number, clobber, debug)
    puts(colored.green(f"---> Successfully imported ({input_vcf}) into {database}"))

@cli.command('merge-batch-vcf', short_help="Combines all sample vcfs databases into a single database")
@click.option('--db-path', '-p', type=click.Path(exists=True), required=True, help="The path to where all the databases for this batch is stored")
@click.option('--cdb', 'caller_db', type=click.Path(), required=True, help="The variant database to import the batch into")
@click.option('--vdb', 'variant_db', type=click.Path(exists=True), required=True, help="The duckdb database to fetch variant ID from")
@click.option('--sdb', 'sample_db', type=click.Path(exists=True), required=True, help="The duckdb database to fetch sample ID from")
@click.option('--caller', 'caller',
              type=click.Choice(['mutect', 'vardict'], case_sensitive=False),
              required=True,
              help="Type of VCF file to import")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def merge_batch_vcf(db_path, caller_db, variant_db, sample_db, caller, batch_number, debug, clobber):
    """
    Ingest the variants in a batch into main variants database
    """
    import ch.vdbtools.importer as importer
    importer.import_caller_batch(db_path, caller_db, variant_db, sample_db, caller, batch_number, debug, clobber)
    log.logit(f"---> Successfully imported variant batch ({batch_number}) into {caller_db}", color="green")

@cli.command('import-sample-variants', short_help="Register the variants for a VCF file into a variant database")
@click.option('--input-vcf', '-i', 'input_vcf', type=click.Path(exists=True), required=True, help="The VCF to be imported")
@click.option('--vdb', 'variant_db', type=click.Path(), required=True, help="The variant database to import the VCF into")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this import set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def import_sample_variants(input_vcf, variant_db, batch_number, debug, clobber):
    """
    Registering variants into the duckdb database.
    """
    import ch.vdbtools.importer as importer
    importer.import_sample_variants(input_vcf, variant_db, batch_number, debug, clobber)
    log.logit(f"---> Successfully imported ({input_vcf})", color="green")

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
    import ch.vdbtools.importer as importer
    importer.import_variant_batch(db_path, variant_db, batch_number, debug, clobber)
    log.logit(f"---> Successfully imported variant batch ({batch_number}) into {variant_db}", color="green")

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
    import ch.vdbtools.dump as dump
    dump.dump_variant_batch(variant_db, header_type, batch_number, chromosome, debug)
    log.logit(f"---> Successfully dumped variant batch ({batch_number}) from {variant_db}", color="green")

# TODO: Write a dump-pileup command to specifically query the pileup and the variants file for variants without any pileup data
@cli.command('dump-variants-pileup', short_help="dumps all variants inside duckdb into a VCF file that needs pileup")
@click.option('--vdb', 'variant_db', type=click.Path(exists=True), required=True, help="The duckdb database to dump the variants from")
@click.option('--pdb', 'pileup_db', type=click.Path(exists=True), required=True, help="The pileup database to check which variants need pileup")
@click.option('--header-type', '-t', type=click.Choice(['simple', 'dummy'], case_sensitive=False), required=True, default="dummy",
                                    help="A pre-existing header type e.g. simple, mutect, vardict, complex, etc.")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--chromosome', '-c', type=click.STRING, default=None, help="The chromosome set of interest")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def dump_variants(variant_db, pileup_db, header_type, batch_number, chromosome, debug):
    """
    Dumps the variants from duckdb into a VCF file
    """
    import ch.vdbtools.dump as dump
    dump.dump_variants_for_pileup(variant_db, pileup_db, header_type, batch_number, chromosome, debug)
    log.logit(f"---> Successfully dumped variant batch ({batch_number}) from {variant_db} that needs pileup", color="green")

@cli.command('import-pon-pileup', short_help="updates variants inside duckdb with PoN pileup information")
@click.option('--vdb', 'variant_db', type=click.Path(exists=True), required=True, help="The duckdb database to fetch variant ID from")
@click.option('--pdb', 'pileup_db', type=click.Path(), required=True, help="The duckdb database to fetch variant ID from")
@click.option('--pon-pileup', '-p', 'pon_pileup', type=click.Path(exists=True), required=True, help="The pon pileup VCF to be imported into the variant database")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def import_pon_pileup(pileup_db, variant_db, pon_pileup, batch_number, debug, clobber):
    """
    Dumps the panel of normal pileup information from into a pileup duckdb
    """
    import ch.vdbtools.importer as importer
    importer.import_pon_pileup(pileup_db, variant_db, pon_pileup, batch_number, debug, clobber)
    log.logit(f"---> Successfully imported PoN Pileup from batch ({batch_number}) into {pileup_db}", color="green")

@cli.command('calculate-fishers-test', short_help="Updates the variants inside Mutect or Vardict tables with p-value from Fisher's Exact Test")
@click.option('--pdb', 'pileup_db', type=click.Path(exists=True), required=True, help="The duckdb database to fetch variant PoN Ref Depth and Alt Depth from")
@click.option('--cdb', 'caller_db', type=click.Path(exists=True), required=True, help="The duckdb database to fetch variant caller information from")
@click.option('--caller', 'caller',
              type=click.Choice(['mutect', 'vardict'], case_sensitive=False),
              required=True,
              help="Type of VCF file to import")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def calculate_fishers_test(pileup_db, caller_db, caller, batch_number, debug):
    """
    Calculates the Fisher's Exact Test for all Variants within the Variant Caller duckdb
    """
    import ch.vdbtools.handlers.callers as callers
    callers.annotate_fisher_test(pileup_db, caller_db, caller, batch_number, debug)
    log.logit(f"---> Successfully calculated the Fisher's Exact Test for variants within ({batch_number}) and {caller_db}", color="green")

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
    import ch.vdbtools.importer as importer
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
    import ch.vdbtools.dump as dump
    dump.dump_variants_for_annotate_pd(annotation_db, batch_number, debug)
    log.logit(f"---> Successfully dumped variant annotations batch ({batch_number}) from {annotation_db}", color="green")

@cli.command('import-annotate-pd', short_help="annotates variants with their pathogenicity")
@click.option('--adb', 'annotation_db', type=click.Path(exists=True), required=True, help="The duckdb database to store the annotation information")
@click.option('--pd', '-p', 'annotate_pd', type=click.Path(exists=True), required=True, help="The CSV File produced from AnnotatePD to be imported into the annotation database")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def import_annotate_pd(annotation_db, annotate_pd, batch_number, debug):
    """
    Dumps the vep information into an annotation duckdb
    """
    import ch.vdbtools.importer as importer
    importer.import_annotate_pd(annotation_db, annotate_pd, batch_number, debug)
    log.logit(f"---> Successfully annotated variants from batch ({batch_number}) in {annotation_db}", color="green")

@cli.command('dump-ch', short_help="Outputs CH Variants from Database")
@click.option('--mcdb', 'mutect_db', type=click.Path(exists=True), required=True, help="The mutect database")
@click.option('--vcdb', 'vardict_db', type=click.Path(exists=True), required=True, help="The vardict database")
@click.option('--adb', 'annotation_db', type=click.Path(exists=True), required=True, help="The annotation database")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def dump_ch_variants(mutect_db, vardict_db, annotation_db, debug):
    """
    Combines all information and outputs CH Variants
    """
    import ch.vdbtools.dump as dump
    dump.dump_ch_variants(mutect_db, vardict_db, annotation_db, debug)
    log.logit(f"---> Successfully dumped CH Variants", color="green")
