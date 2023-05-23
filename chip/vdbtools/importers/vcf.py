import vcfpy
import collections
import pysam
import os, io, gzip

import importlib.resources
import chip.utils.logger as log
import pandas as pd

from clint.textui import indent, puts_err, puts

class Vcf:
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
        self.reader = vcfpy.Reader.from_path(self.vcf_path)

    def reset_reader(self):
        self.reader = vcfpy.Reader.from_path(self.vcf_path)

def duckdb_load_df_file(duckdb_connection, df, table):
    sql = f"""
        INSERT INTO {table} SELECT * FROM df
    """
    log.logit(f"Starting to insert pandas dataframe into duckdb")
    duckdb_connection.execute(sql)
    log.logit(f"Finished inserting pandas dataframe into duckdb")

def vcf_to_pd(input_vcf, what_process, batch_number, debug):
    log.logit(f"Processing: {input_vcf}")
    dispatch = {
        'variants'  : variants_to_df,
        'caller'    : caller_to_df,
        'pileup'    : pileup_to_df
    }
    function = dispatch[what_process]
    return function(input_vcf, batch_number, debug)

def variants_to_df(input_vcf, batch_number, debug):
    log.logit(f"Reading in the Variants VCF...")
    res = pd.read_csv(input_vcf,
            comment='#',
            compression='gzip',
            sep='\t',
            header=None,
            usecols=[0,1,3,4]).rename(columns={0: "chrom",
                                              1: "pos",
                                              3: "ref",
                                              4: "alt"})
    log.logit(f"Finished reading in the variant VCF")
    res['snp'] = (res['ref'].apply(len) == 1) & (res['alt'].apply(len) == 1)
    res['qc_pass'] = None
    res['batch'] = batch_number
    res['start'] = res['pos']
    res['end'] = res['pos'] + res['alt'].apply(len)
    res['key'] = res['chrom'] + ':' + res['pos'].astype(str) + ':' + res['ref'] + ':' + res['alt']
    total = len(res)
    return total, res

def pileup_to_df(input_vcf, batch_number, debug):
    log.logit(f"Reading in the Pileup VCF...")
    res = pd.read_csv(input_vcf,
            comment='#',
            compression='gzip',
            sep='\t',
            header=None,
            usecols=[0,1,3,4,7]).rename(columns={0: "chrom",
                                              1: "pos",
                                              3: "ref",
                                              4: "alt",
                                              7: "info"})
    log.logit(f"Finished reading in the pileup VCF")
    log.logit(f"Formatting the dataframe...")
    res['key'] = res['chrom'] + ':' + res['pos'].astype(str) + ':' + res['ref'] + ':' + res['alt']
    res[['PoN_RefDepth','PoN_AltDepth']] = res['info'].str.split(";", expand=True)
    res['PoN_RefDepth'] = res['PoN_RefDepth'].str.split("=").str[1]
    res['PoN_AltDepth'] = res['PoN_AltDepth'].str.split("=").str[1]
    res['batch'] = batch_number
    res['variant_id'] = None
    res = res[['key', 'PoN_RefDepth', 'PoN_AltDepth', 'batch', 'variant_id']]
    log.logit(f"Finished preparing the dataframe...")
    total = len(res)
    return total, res

def getFields(fields):
    typeMap = {
        'String': str,
        'Float': float,
        'Integer': int,
        'Flag': bool
    }
    fields = [[x.split('=')[-1] for x in field[:3]] for field in fields]
    for field in fields:
        if field[1] == 'R':
            field[2] = str
        elif field[2] == 'Integer' and field[1] != 'A' and field[1] != 'G' and int(field[1]) > 1:
            field[2] = str
        else:
            field[2] = typeMap[field[2]]
    fields = {field[0]:field[2] for field in fields}
    return fields

# def splitInfor(x, columnType):
#     res = {k: columnType[k](v) for k, v in [y.split('=') if '=' in y else [y, True] for y in x.split(';')]}
#     return [res[col] if col in res else None for col in columnType]

def splitInfor(info_tags, fields):
    x = [tag.split('=') if '=' in tag else [tag, True] for tag in info_tags.split(';')] # Splits the Tags into Key = Value or Flag = True
    res = {k: fields[k](v) for k, v in x}                                               # Converts the Value into the correct Type (via Fields Mapping)
    return [res[field] if field in res else None for field in fields]                   # For each "info_tag/key" (that should be here), assign it the value to the column or else None

def combineAndSplit(format, sample, fields):
    res = {format.split(':')[field]: sample.split(':')[field] for field in range(len(format.split(':')))}   # For each Format Flag, get the dictionary pairing e.g. GT: '0/1'
    return [fields[field](res[field]) if field in res else None for field in fields]                        # For each "format_tag/key", assign it the value after mapping it to the correct type

def caller_to_df(input_vcf, batch_number, debug):
    log.logit(f"Reading in the VCF...")
    info_fields, format_fields = [], []
    with gzip.open(input_vcf, 'rt') as f:
        for line in f:
            if line.startswith('##INFO'):
                info_fields.append(line.split(','))
            elif line.startswith('##FORMAT'):
                format_fields.append(line.split(','))
            elif line.startswith('#CHROM'):
                header = line.strip('\n').split('\t')
                break
    f.close()
    info_fields = getFields(info_fields)
    format_fields = getFields(format_fields)
    header[-1] = "SAMPLE"

    res = pd.read_csv(input_vcf,
                comment='#',
                compression='gzip',
                sep='\t',
                header=None,
                names=header).rename(columns={'#CHROM':'CHROM'})
    log.logit(f"Finished reading in the VCF")
    log.logit(f"Parsing and formatting INFO and FORMAT columns...")
    info_field_values = pd.DataFrame(list(res['INFO'].apply(lambda x: splitInfor(x, info_fields) ).values),
                                    columns=info_fields.keys())
    info_field_values.columns = info_field_values.columns.str.lower()
    info_field_values = info_field_values.add_prefix('info_')
    format_field_values = pd.DataFrame(list(res[['FORMAT', 'SAMPLE']].apply(lambda x: combineAndSplit(*x, format_fields), axis=1).values),
                                    columns=format_fields.keys())
    format_field_values.columns = format_field_values.columns.str.lower()
    format_field_values = format_field_values.add_prefix('format_')

    res = pd.concat([res.drop('INFO', axis=1), info_field_values], axis=1)
    res = pd.concat([res.drop(['FORMAT', 'SAMPLE'], axis=1), format_field_values], axis=1)
    log.logit(f"Finished formatting INFO and FORMAT columns...")
    res['key'] = res['CHROM'] + ':' + res['POS'].astype(str) + ':' + res['REF'] + ':' + res['ALT'] #Key needed to get VariantID
    res['FILTER'] = res['FILTER'].str.split(";")
    total = len(res)
    return total, res

def load_simple_header(header_type):
    # TODO
    if header_type == "complex":
        log.logit("Loading complex VCF header...")
        return(source)
    elif header_type == "simple":
        log.logit("Loading simple VCF header...")
        source = importlib.resources.files('chip.resources.vcf').joinpath('simple.header')
        return(source)
    else:
        log.logit("Loading default VCF header...")
        source = importlib.resources.files('chip.resources.vcf').joinpath('dummy.header')
        return(source)

def variants_to_vcf(duckdb_variants, header_type, batch_number, chromosome, debug):
    if chromosome == None:
        filename = f"batch-" + str(batch_number) + ".vcf"
    else:
        filename = chromosome + '.vcf'
    reader = open(load_simple_header(header_type), "r")
    header = reader.read()
    reader.close()
    log.logit(f"Preparing to write the variants to the VCF: {filename}")
    log.logit(f"Starting to sort variants...")
    df = duckdb_variants.df().sort_values(by=['chrom', 'pos'], ascending = True)
    df['qual'] = "."
    df['filter'] = "PASS"
    df['info'] = "."
    df['format'] = "GT"
    df['sample'] = "0/1"
    df = df[['chrom', 'pos', 'variant_id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'sample']]
    total = len(df)
    with open(filename, 'w') as writer:
        writer.write(header)
    log.logit(f"Writing to VCF...")
    df.to_csv(filename, sep="\t", mode='a', index=False, header=False)
    log.logit(f"Finished writing {total} variants to the VCF: {filename}")
    log.logit(f"Compressing and indexing the VCF: {filename}")
    pysam.tabix_compress(filename, filename + '.gz', force = True)
    os.unlink(filename)
    pysam.tabix_index(filename + '.gz', preset="vcf", force = True)
