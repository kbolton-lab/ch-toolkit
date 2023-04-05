import os, csv
import vcfpy

import chip.utils.logger as log
from clint.textui import indent, puts_err, puts

def in_lsf_session():
    result = True if 'LSB_JOBID' in os.environ else False
    return result

def create_tmp_csv(work_dir, batch_number, chromosome):
    log.logit(f"Creating temporary CSV file in: {work_dir}")
    tmpdir = work_dir
    filename = None
    if in_lsf_session() and tmpdir == '/tmp':
        tmpdir = os.path.join('/tmp', f"{os.environ['LSB_JOBID']}.tmpdir")
    if chromosome is not None:
        filename = f"batch-{batch_number}.{chromosome}.csv"
    else:
        filename = f"batch-{batch_number}.csv"

    tmp_path = os.path.join(tmpdir, filename)
    if os.path.exists(tmp_path):
        log.logit(f"Deleting existing temporary csv file: {tmp_path}")
        os.unlink(tmp_path)

    return tmp_path

def duckdb_load_csv_file(duckdb_connection, temp_csv, table):
    sql = f"""
        COPY {table} FROM '{temp_csv}' (AUTO_DETECT TRUE)
    """
    log.logit(f"Starting to load csv into duckdb")
    duckdb_connection.execute(sql)
    log.logit(f"Finished loading csv into duckdb")

def vcf_to_csv(input_vcf, what_process, headers, batch_number, work_dir, debug):
    log.logit(f"Processing: {input_vcf}")
    dispatch = {
        'variants'  : variants_to_csv,
        'caller'    : caller_to_csv
    }
    function = dispatch[what_process]
    return function(input_vcf, headers, batch_number, work_dir, debug)

def variants_to_csv(input_vcf, headers, batch_number, work_dir, debug):
    reader = vcfpy.Reader.from_path(input_vcf)
    res = []
    with indent(4, quote=' >'):
        for record in reader:
            (chrom, pos, ref, alt) = (record.CHROM, record.POS, record.REF, record.ALT[0].value)
            start = pos
            stop = start + len(alt)
            snp = True if len(ref) == len(alt) == 1 else False
            item = {
                'chrom'      : chrom,
                'pos'        : pos,
                'ref'        : ref,
                'alt'        : alt,
                'snp'        : snp,
                'qc_pass'    : None,
                'batch'      : batch_number,
                'start'      : start,
                'stop'       : stop,
                'PoN_RefDepth': None,
                'PoN_AltDepth': None
            }
            if debug: log.logit(f"{item}")
            res.append(item)
    count = len(res)
    tmp_csv = create_tmp_csv(work_dir, batch_number, None)
    log.logit(f"Converting variants in VCF file to CSV file: {tmp_csv}")
    with open(tmp_csv, 'wt') as fz:
        writer = csv.DictWriter(fz, fieldnames=headers)
        writer.writeheader()
        writer.writerows(res)
    fz.close()
    log.logit(f"Finished converting variants from: {input_vcf} into csv file. {count} total variants")
    return(count, tmp_csv)

def caller_to_csv():
    return None
