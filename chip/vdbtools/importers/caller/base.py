import abc, sqlite3

import chip.vdbtools.importers.vcf as vcf

class BaseIngestor(metaclass=abc.ABCMeta):
    def __init__(self, db_path, input_vcf):
        self.db_path = db_path
        self.db_connection = sqlite3.connect(self.db_path)
        self.vcf = vcf.Vcf(input_vcf)

    @abc.abstractmethod
    def ingest_vcf(self):
        pass

    @abc.abstractmethod
    def resolve_variant(self, record):
        pass
