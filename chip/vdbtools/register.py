import sys, time

import chip.utils.logger as log

from clint.textui import indent, puts_err, puts
import redis
import vcfpy

def _redis_connect(host, port):
    r = redis.Redis(host=host, port=port, db=0)
    return r

def _process_vcf(input_vcf, redis_db, batch_number, debug):
    puts(f"Processing: {input_vcf}")
    reader = vcfpy.Reader.from_path(input_vcf)
    (total, new, seen) = (0, 0, 0)
    redis_set = f"batch:{batch_number}"
    for record in reader:
        if total % 1000 == 0:
            log.logit(f"variants processed - total: {total} | new: {new} | seen: {seen}")
        (chrom, pos, ref, alt) = (record.CHROM, str(record.POS), record.REF, record.ALT[0].value)
        key = ':'.join([chrom, pos, ref, alt])
        if redis_db.exists(key):
            if debug: puts_err(f"variant: '{key}' already seen -- skipping")
            seen += 1
        else:
            time.sleep(2)
            pipe = redis_db.pipeline()
            pipe.sadd(redis_set, key)
            pipe.incr('variant_id')
            pipe.copy('variant_id', key)
            vals = pipe.execute()
            variant_id = redis_db.get(key)
            if debug: puts_err(f"variant: '{key}' => {variant_id} [{vals}]")
            new += 1
        total += 1
    return { 'total' : total, 'new' : new, 'seen' : seen }

def check_global_variant_counter(redis_db, start=0):
    if not redis_db.exists('variant_id'):
        log.logit(f"Creating the variant id global counter: variant_id={start}", color="yellow")
        redis_db.set('variant_id', start)

def import_vcf(input_vcf, redis_host, redis_port, batch_number, debug):
    log.logit(f"Registering variants from file: {input_vcf}", color="green")
    with indent(4, quote=' >'):
        redis_db = _redis_connect(redis_host, redis_port)
        check_global_variant_counter(redis_db)
        counts = _process_vcf(input_vcf, redis_db, batch_number, debug)
    log.logit(f"Finished registering variants")
    log.logit(f"Variants Processed - Total: {counts['total']} | New: {counts['new']} | Seen: {counts['seen']}", color="green")
    log.logit(f"All Done!", color="green")
