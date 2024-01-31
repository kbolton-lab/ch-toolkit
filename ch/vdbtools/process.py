import sys
import multiprocessing as mp
import ch.utils.logger as log

def annotate_fisher_test(pileup_db, caller_db, caller, batch_number, by_chromosome, debug):
    import ch.vdbtools.handlers.callers as callers
    callers.annotate_fisher_test(pileup_db, caller_db, caller, batch_number, by_chromosome, debug)

def db_to_chromosome(db, which_db, batch_number, chromosome, cores, debug):
    dispatch = {
        'mutect'  : caller_to_chromosome,
        'vardict' : caller_to_chromosome,
        'variant' : variant_to_chromosome,
        'annotation' : annotation_to_chromosome
    }
    function = dispatch[which_db]
    function(db, which_db, batch_number, chromosome, cores, debug)

def set_options(db, batch_number, chromosome, debug):
    if batch_number is None:
        if chromosome is None:
            log.logit(f"Splitting {db} variants into individual chromosomes", color="green")
            chromosome = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 
               'chr21', 'chr22', 'chrX', 'chrY']
        else:
            log.logit(f"Splitting {db} variants into chromosome: {chromosome}", color="green")
            chromosome = [chromosome]
        base_db = db.replace('.db', '')
    else:
        if chromosome is None:
            log.logit(f"Splitting {db} variants into individual chromosomes for batch: {batch_number}", color="green")
            chromosome = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 
               'chr21', 'chr22', 'chrX', 'chrY']
        else:
            log.logit(f"Splitting {db} variants into chromosome: {chromosome} for batch: {batch_number}", color="green")
            chromosome = [chromosome]
        base_db = db.replace('.db', '') + f".batch{batch_number}"
    return batch_number, chromosome, base_db

def caller_to_chromosome(caller_db, caller, batch_number, chromosome, cores, debug):
    batch_number, chromosome, base_db = set_options(caller_db, batch_number, chromosome, debug)
    import ch.vdbtools.handlers.callers as callers
    with mp.Pool(cores) as p:
        p.starmap(callers.caller_to_chromosome, [(caller_db, caller, batch_number, chrom, base_db, debug) for chrom in chromosome])
    #callers.caller_to_chromosome(caller_db, caller, batch_number, chromosome, base_db, debug)

def variant_to_chromosome(variant_db, variant, batch_number, chromosome, cores, debug):
    batch_number, chromosome, base_db = set_options(variant_db, batch_number, chromosome, debug)
    import ch.vdbtools.handlers.variants as variants
    with mp.Pool(cores) as p:
        p.starmap(variants.variant_to_chromosome, [(variant_db, variant, chrom, base_db, debug) for chrom in chromosome])
    #variants.variant_to_chromosome(variant_db, variant, chromosome, base_db, debug)

def annotation_to_chromosome(annotation_db, annotation, batch_number, chromosome, cores, debug):
    batch_number, chromosome, base_db = set_options(annotation_db, batch_number, chromosome, debug)
    import ch.vdbtools.handlers.annotations as annotations
    with mp.Pool(cores) as p:
        p.starmap(annotations.annotation_to_chromosome, [(annotation_db, annotation, chrom, base_db, debug) for chrom in chromosome])
    #annotations.annotation_to_chromosome(annotation_db, annotation, chromosome, base_db, debug)

def pileup_to_chromosome(pileup_db, pileup, batch_number, chromosome, cores, debug):
    batch_number, chromosome, base_db = set_options(pileup_db, batch_number, chromosome, debug)
    import ch.vdbtools.handlers.variants as variants
    with mp.Pool(cores) as p:
        p.starmap(variants.pileup_to_chromosome, [(pileup_db, pileup, chrom, base_db, debug) for chrom in chromosome])
    #variants.pileup_to_chromosome(pileup_db, pileup, chromosome, base_db, debug)

def chromosome_to_caller(chr_path, caller_db, caller, debug):
    import ch.vdbtools.handlers.callers as callers
    callers.chromosome_to_caller(chr_path, caller_db, caller, debug)

def ch_variants_only(caller_db, caller, annotation_db, batch_number, chromosome, cores, debug):
    batch_number, chromosome, base_db = set_options(caller_db, batch_number, chromosome, debug)
    import ch.vdbtools.analysis.ch as ch
    #for chrom in chromosome:
    #    ch.ch_variants_only(mutect_db, vardict_db, annotation_db, chrom, debug)
    with mp.Pool(cores) as p:
        p.starmap(ch.ch_variants_only, [(caller_db, base_db, caller, annotation_db, chrom, debug) for chrom in chromosome])
    ch.merge_ch_variants(base_db, caller)