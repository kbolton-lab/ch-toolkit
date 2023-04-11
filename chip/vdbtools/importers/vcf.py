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
    res['PoN_RefDepth'] = None
    res['PoN_AltDepth'] = None
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
    res[['PoN_RefDepth','PoN_AltDepth']] = res['info'].str.split(";", expand=True)
    res['PoN_RefDepth'] = res['PoN_RefDepth'].str.split("=").str[1]
    res['PoN_AltDepth'] = res['PoN_AltDepth'].str.split("=").str[1]
    res['snp'] = None
    res['qc_pass'] = None
    res['batch'] = batch_number
    res['start'] = None
    res['end'] = None
    res['key'] = res['chrom'] + ':' + res['pos'].astype(str) + ':' + res['ref'] + ':' + res['alt']
    res = res[['chrom', 'pos', 'ref', 'alt', 'snp', 'qc_pass', 'batch', 'start', 'end', 'PoN_RefDepth', 'PoN_AltDepth', 'key']]
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

# Input: info_tags row     [['AS_FilterStatus', 'weak_evidence,map_qual']
#                          ['AS_SB_TABLE', '7,13|1,2'] ... etc ]
# Output: list of values   ['weak_evidence,map_qual', '7,13|1,2', None, None, ... etc ]
def splitInfor(x, fields):
    x = [y.split('=') for y in x.split(';')]
    res = {}
    for kv in x:
        if len(kv) == 2:
            res[kv[0]] = kv[1]
        else:
            res[kv[0]] = True
    res = [fields[field](res[field]) if field in res else None for field in fields]
    return res

def combineAndSplit(format, sample, fields):
    res = {format.split(':')[field]: sample.split(':')[field] for field in range(len(format.split(':')))}
    res = [fields[field](res[field]) if field in res else None for field in fields]
    return res

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
    variants = duckdb_variants.df().sort_values(by=['chrom', 'pos'], ascending = True)
    variants['ROWID']
    variants['qual'] = "."
    variants['filter'] = "PASS"
    variants['info'] = "."
    variants['format'] = "GT"
    variants['sample'] = "0/1"
    variants = variants[['chrom', 'pos', 'ROWID', 'ref', 'alt', 'qual', 'filter', 'info', 'format']]
    total = len(variants.index)
    with open(filename, 'w') as writer:
        writer.write(header)
    log.logit(f"Writing to VCF...")
    variants.to_csv(filename, sep="\t", mode='a', index=False, header=False)
    log.logit(f"Finished writing {total} variants to the VCF: {filename}")
    log.logit(f"Compressing and indexing the VCF: {filename}")
    pysam.tabix_compress(filename, filename + '.gz', force = True)
    os.unlink(filename)
    pysam.tabix_index(filename + '.gz', preset="vcf", force = True)

# DEPRECATED
def write_variants_to_vcf(duckdb_variants, header_type, batch_number, chromosome, debug):
    reader = vcfpy.Reader.from_path(load_simple_header(header_type))
    header = reader.header
    #header = vcfpy.Header(samples = str(batch_number))
    #meta_line = vcfpy.MetaHeaderLine('fileformat', '', {'META': 'VCFv4.3'})
    #header.add_line(meta_line)
    # Add the fileformat meta header line
    #header.add_line(vcfpy.MetaHeaderLine('fileformat', None, {'META': 'VCFv4.3'}))

    #info_field = vcfpy.OrderedDict([('ID', 'FAKE'), ('Number', '1'), ('Type', 'String'), ('Description', 'Dummy Value')])
    #header.add_info_line(vcfpy.OrderedDict(info_field))
    #format_field = vcfpy.OrderedDict([('ID', 'GT'), ('Number', '1'), ('Type', 'String'), ('Description', 'Genotype')])
    #header.add_format_line(vcfpy.OrderedDict(format_field))
    header.samples = vcfpy.SamplesInfos(str(batch_number))
    if chromosome == None:
        filename = f"batch-" + str(batch_number) + ".vcf"
    else:
        filename = chromosome + '.vcf'
    writer = vcfpy.Writer.from_path(filename, header)
    log.logit(f"Preparing to write the variants to the VCF: {filename}")
    with indent(4, quote=' >'):
        log.logit(f"Starting to sort variants...")
        variants = sorted(duckdb_variants.fetchall(), key=lambda row: (row[1], row[2]))
        total = len(variants)
        log.logit(f"Finished sorting {total} variants")
        log.logit(f"Writing to VCF...")
        for i, variant in enumerate(variants):
            record = vcfpy.Record(
                CHROM = variant[1],
                POS = variant[2],
                ID = [str(variant[0])],
                REF = variant[3],
                ALT = [vcfpy.Substitution(type_ = "", value = variant[4])],
                QUAL = ".",
                FILTER = ["PASS"],
                INFO = {"FAKE" : "."},
                FORMAT = ["GT"],
                calls = [vcfpy.Call(sample = str(batch_number), data = {"GT" : "0/1"})]
            )
            if debug: log.logit(str(record))
            writer.write_record(record)
            if i % 1_000_000 == 0 and i != 0:
                log.logit(f"{i} variants written to: {filename}")
    log.logit(f"Finished writing {total} variants to the VCF: {filename}")
    writer.close()
    log.logit(f"Compressing and indexing the VCF: {filename}")
    pysam.tabix_compress(filename, filename + '.gz', force = True)
    os.unlink(filename)
    pysam.tabix_index(filename + '.gz', preset="vcf", force = True)
