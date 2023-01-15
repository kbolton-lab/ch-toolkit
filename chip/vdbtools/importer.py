import os, path

def import_vcf(db_path, input_vcf, caller):
    dispatch = {
        'lofreq'  : import_lofreq,
        'mutect'  : import_mutect,
        'vardict' : import_vardict,
        'pindel'  : import_pindel,
    }
    pass

def import_lofreq():
    pass

def import_mutect():
    pass

def import_vardict():
    pass

def import_pindel():
    pass
