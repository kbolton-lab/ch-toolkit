import os, sys, signal
import click
from clint.textui import puts, colored
import ch.utils.database as db

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
@click.option('--threads', 'cores', type=click.INT, required=False, show_default=True, default=1, help="Number of Threads used for parallelization")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
@click.option('--clobber', '-f', is_flag=True, show_default=True, default=False, required=True, help="If exists, delete existing duckdb file and then start from scratch")
def merge_batch_vcf(db_path, caller_db, variant_db, sample_db, caller, batch_number, cores, debug, clobber):
    """
    Ingest the variants in a batch into main variants database
    """
    import ch.vdbtools.importer as importer
    importer.import_caller_batch(db_path, caller_db, variant_db, sample_db, caller, batch_number, cores, debug, clobber)
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

@cli.command('bcbio-filter', short_help="Filters variants in the Vardict Database using the BCBIO filter")
@click.option('--vcdb', 'vardict_db', type=click.Path(), required=True, help="Vardict database to calculate the BCBIO filter")
@click.option('--recalculate', '-r', is_flag=True, show_default=True, default=False, required=True, help="Recalculates the BCBIO Filter String Automatically (Ignores Default Parameters)")
@click.option('--low-dp-af', 'low_depth_for_allele_frequency', type=click.INT, required=True, show_default=True, default=6, help="Filter cut-off for regions with low coverage for allele frequency")
@click.option('--total-depth', '--dp', type=click.INT, required=True, show_default=True, default=10, help="Cut-off for what depth is considered low coverage")
@click.option('--mean-quality-score', '--qual', type=click.INT, required=True, show_default=True, default=30, help="Cut-off for low quality region scores")
@click.option('--batch-number', '-b', type=click.INT, required=True, help="The batch number of this variant set")
@click.option('--by_chromosome', '-c', is_flag=True, show_default=True, default=False, required=True, help="By chromosome or all at once")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def bcbio_filter(vardict_db, recalculate, low_depth_for_allele_frequency, total_depth, mean_quality_score, batch_number, by_chromosome, debug):
    """
    Performs the BCBIO Filter on the Vardict Database: https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/variation/vardict.py#L251\n\n
    The variant_calling.wdl should automatically perform the BCBIO filter. However, if the filter_string was not properly set, this filter can be run at this step.\n
    As a rule of thumb, parameters for low-dp-af, total-depth, and mean-quality-score are calculated using a small subset then identifying the cut-off for 2% of the left side samples.\n
    However, if this was not done previously, the --recalculate parameter can be used to do this automatically.\n\n
    If the user knows what values to use set them using the follow options.
    """
    vardict_connection = db.duckdb_connect_ro(vardict_db)
    check_vardict = vardict_connection.execute(f"SELECT * FROM information_schema.tables WHERE table_name = 'vardict'").fetchall()
    vardict_connection.close()
    if len(check_vardict) > 0:
        import ch.vdbtools.process as process
        if recalculate:
            low_depth_for_allele_frequency, total_depth, mean_quality_score = process.recalculate_bcbio_parameters(vardict_db, low_depth_for_allele_frequency, debug)
        log.logit(f"The BCBIO Filter String is: ((FMT/AF * FMT/DP < {low_depth_for_allele_frequency}) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 2.0) || (FMT/DP < {total_depth}) || (INFO/QUAL < {mean_quality_score})))")
        process.bcbio_filter(vardict_db, low_depth_for_allele_frequency, total_depth, mean_quality_score, batch_number, by_chromosome, debug)
        log.logit(f"---> Successfully performed BCBIO filtering of regions with low coverage for allele fractions within (batch {batch_number}) in {vardict_db}", color="green")
    else:
        log.logit("ERROR: BCBIO Filter is only needed for the Vardict Database. Please provide the vardict database for --vcdb", color="red")

#! Convert a database into chromosomes
@cli.command('database-to-chromosome', short_help="Splits <Mutect|Vardict|Variant|Annotation|Pileup> database into individual chromosomes")
@click.option('--db', 'db', type=click.Path(), required=True, help="The database to split into chromosomes")
@click.option('--which_db', 'which_db',
              type=click.Choice(['mutect', 'vardict', 'variant', 'annotation', 'pileup'], case_sensitive=False),
              required=True,
              help="The specific database being processed")
@click.option('--batch-number', '-b', type=click.INT, default=None, help="The batch number in case only want chromosome database for a subset")
@click.option('--threads', 'cores', type=click.INT, required=False, show_default=True, default=1, help="Number of Threads used for parallelization")
@click.option('--chromosome', '-c', type=click.STRING, default=None, help="The chromosome set of interest")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def caller_to_chromosome(db, which_db, batch_number, chromosome, cores, debug):
    """
    Splits the <Mutect|Vardit|Variant|Annotation|Pileup> database into individual chromosomes\n
    Has the functionality to be specific about which batch needs to be processed as well as which chromosome\n
    Default is to split the database into individual chromosomes (1-22, X, Y) for all batches
    variant_db, pileup_db, and annotation_db CANNOT be split by batch - ONLY by chromosome
    """
    import ch.vdbtools.process as process
    process.db_to_chromosome(db, which_db, batch_number, chromosome, cores, debug)
    if batch_number != None:
        if chromosome != None:
            log.logit(f"---> Successfully split {db} from batch: {batch_number} into chromosome: {chromosome}", color="green")
        else:
            log.logit(f"---> Successfully split {db} from batch: {batch_number} into ALL chromosomes", color="green")
    else:
        if chromosome != None:
            log.logit(f"---> Successfully split {db} into chromosome: {chromosome}", color="green")
        else:
            log.logit(f"---> Successfully split {db} into ALL chromosomes", color="green")

# TODO: Convert chromosome split database back into a single caller database
@cli.command('chromosome-to-caller', short_help="Combines all chromosome databases into a single <Mutect|Vardict> database")
@click.option('--chr-path', '-p', type=click.Path(exists=True), required=True, help="The path to where all the databases to be combined are stored")
@click.option('--cdb', 'caller_db', type=click.Path(), required=True, help="The caller database to combine the chromosomes")
@click.option('--caller', 'caller',
              type=click.Choice(['mutect', 'vardict'], case_sensitive=False),
              required=True,
              help="Type of VCF file to import")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def chromosome_to_caller(chr_path, caller_db, caller, debug):
    """
    After splitting the caller databases into individual variants to do processing, combine the information into a central system\n
    NOTE: This will create a completely new <Mutect|Vardict> database. If columns have been changed in the chromosomes, then the new database will reflect those changes.
    """
    import ch.vdbtools.process as process
    process.caller_to_chromosome(chr_path, caller_db, caller, debug)
    log.logit(f"---> Successfully imported chromosomes from {chr_path} into {caller_db}", color="green")

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
@click.option('--by_chromosome', '-c', is_flag=True, show_default=True, default=False, required=True, help="By chromosome or all at once")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def calculate_fishers_test(pileup_db, caller_db, caller, batch_number, by_chromosome, debug):
    """
    Calculates the Fisher's Exact Test for all Variants within the Variant Caller
    """
    import ch.vdbtools.process as process
    process.annotate_fisher_test(pileup_db, caller_db, caller, batch_number, by_chromosome, debug)
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
@click.option('--prefix', '-p', type=click.STRING, default="ch_pd", help="The output prefix e.g. <prefix>.all.csv")
@click.option('--pvalue', '-v', type=click.FLOAT, default=1.260958e-09, help="The p-value cut-off value for the Fisher's exact test for the PoN")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def dump_ch_variants(mutect_db, vardict_db, annotation_db, prefix, pvalue, debug):
    """
    Combines all information and outputs CH Variants
    """
    import ch.vdbtools.dump as dump
    dump.dump_ch_variants(mutect_db, vardict_db, annotation_db, prefix, pvalue, debug)
    log.logit(f"---> Successfully dumped CH Variants", color="green")

@cli.command('reduce-db', short_help="Reduces the size of the mutect_db and vardict_db databases to only CH possible variants")
@click.option('--cdb', 'caller_db', type=click.Path(exists=True), required=True, help="The mutect or vardict database")
@click.option('--caller', 'caller',
              type=click.Choice(['mutect', 'vardict'], case_sensitive=False),
              required=True,
              help="Select between: mutect or vardict")
@click.option('--adb', 'annotation_db', type=click.Path(exists=True), required=True, help="The annotation database")
@click.option('--batch-number', '-b', type=click.INT, default=None, help="The batch number in case only want chromosome database for a subset")
@click.option('--threads', 'cores', type=click.INT, required=False, show_default=True, default=1, help="Number of Threads used for parallelization")
@click.option('--chromosome', '-c', type=click.STRING, default=None, help="The chromosome set of interest")
@click.option('--debug', '-d', is_flag=True, show_default=True, default=False, required=True, help="Print extra debugging output")
def dump_ch_variants(caller_db, caller, annotation_db, batch_number, chromosome, cores, debug):
    """
    Reduces the size of mutect_db and vardict_db to only include possible CH variants
    """
    import ch.vdbtools.process as process
    process.ch_variants_only(caller_db, caller, annotation_db, batch_number, chromosome, cores, debug)
    if batch_number != None:
        if chromosome != None:
            log.logit(f"---> Successfully reduced {caller_db} from batch: {batch_number} for chromosome: {chromosome}", color="green")
        else:
            log.logit(f"---> Successfully reduced {caller_db} from batch: {batch_number} for ALL chromosomes", color="green")
    else:
        if chromosome != None:
            log.logit(f"---> Successfully reduced {caller_db} for chromosome: {chromosome}", color="green")
        else:
            log.logit(f"---> Successfully reducced {caller_db} for ALL chromosomes", color="green")
