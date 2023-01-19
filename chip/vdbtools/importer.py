import sys

def import_vcf(db_path, input_vcf, caller):
    dispatch = {
        'lofreq'  : import_lofreq,
        'mutect'  : import_mutect,
        'vardict' : import_vardict,
        'pindel'  : import_pindel,
    }
    function = dispatch[caller]
    function(db_path, input_vcf)

def import_lofreq(db_path, input_vcf):
    import chip.vdbtools.importers.caller.lofreq as lofreq
    ingestor = lofreq.Ingestor(db_path, input_vcf)
    ingestor.ingest_vcf()

def import_mutect(db_path, input_vcf):
    sys.exit("[err] Please implement me -- import_mutect !")

def import_vardict(db_path, input_vcf):
    sys.exit("[err] Please implement me -- import_vardict !")

def import_pindel(db_path, input_vcf):
    sys.exit("[err] Please implement me -- import_pindel !")
