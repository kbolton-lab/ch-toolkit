import vcfpy

class Vcf:
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
        self.reader = vcfpy.Reader.from_path(self.vcf_path)

    def reset_reader(self):
        self.reader = vcfpy.Reader.from_path(self.vcf_path)
