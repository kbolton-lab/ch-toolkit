import sqlite3

import chip.vdbtools.importers.caller.base as base

class Ingestor(base.BaseIngestor):
    def __init__(self, db_path, input_vcf):
        super().__init__(db_path, input_vcf)

    def ingest_vcf(self):
        for record in self.vcf.reader:
            variant = self.resolve_variant(record)

    def resolve_variant(self, record):
        pass
