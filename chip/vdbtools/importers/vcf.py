import vcfpy
import collections
import pysam
import os
import importlib.resources
from clint.textui import indent, puts_err, puts
import chip.utils.logger as log

class Vcf:
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
        self.reader = vcfpy.Reader.from_path(self.vcf_path)

    def reset_reader(self):
        self.reader = vcfpy.Reader.from_path(self.vcf_path)

def load_simple_header(header_type):
    # TODO
    if header_type == "complex":
        log.logit("Loading complex VCF header...")
        return(source)
    else:
        log.logit("Loading default VCF header...")
        source = importlib.resources.files('chip.resources.vcf').joinpath('dummy.header')
        return(source)

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
        filename = f"batch_" + str(batch_number) + ".vcf"
    else:
        filename = chromosome + '.vcf'
    writer = vcfpy.Writer.from_path(filename, header)
    log.logit(f"Preparing to write the variants to the VCF: {filename}")
    with indent(4, quote=' >'):
        #>> [x for x in a if x <= 5]
        for variant in sorted(duckdb_variants.fetchall(), key=lambda row: (row[1], row[2])):
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
    log.logit(f"Finished writing the variants to the VCF: {filename}")
    writer.close()
    log.logit(f"Compressing and indexing the VCF: {filename}")
    pysam.tabix_compress(filename, filename + '.gz', force = True)
    os.unlink(filename)
    pysam.tabix_index(filename + '.gz', preset="vcf", force = True)
