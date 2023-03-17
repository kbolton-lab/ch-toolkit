import vcfpy
import collections
import pysam

class Vcf:
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
        self.reader = vcfpy.Reader.from_path(self.vcf_path)

    def reset_reader(self):
        self.reader = vcfpy.Reader.from_path(self.vcf_path)

def write_variants_to_vcf(duckdb_variants, header, batch_number, chromosome):
    reader = vcfpy.Reader.from_path(header)
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
    writer = vcfpy.Writer.from_path(chromosome + '.vcf', header)
    for variant in sorted(duckdb_variants.fetchall(), key=lambda row: (row[0], row[1])):
        record = vcfpy.Record(
            CHROM = variant[0],
            POS = variant[1],
            ID = ".",
            REF = variant[2],
            ALT = [vcfpy.Substitution(type_ = "", value = variant[3])],
            QUAL = ".",
            FILTER = ["PASS"],
            INFO = {"FAKE" : "."},
            FORMAT = ["GT"],
            calls = [vcfpy.Call(sample = str(batch_number), data = {"GT" : "0/1"})]
        )
        writer.write_record(record)
    writer.close()
    pysam.tabix_compress(chromosome + '.vcf', chromosome + '.vcf.gz', force = True)
    pysam.tabix_index(chromosome + '.vcf.gz', preset="vcf", force = True)
