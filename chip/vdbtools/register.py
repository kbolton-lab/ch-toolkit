import sys

from clint.textui import puts, indent
import redis
import vcfpy

def _redis_connect(host, port):
    r = redis.Redis(host=host, port=port, db=0)
    return r

def _process_vcf(input_vcf, redis_db, batch_number):
    puts(f"Processing: {input_vcf}")
    reader = vcfpy.Reader.from_path(vcf_file)
    counter = 1
    redis_set = f"new:{batch_number}"
    for record in reader:
        if counter % 100 == 0:
            puts(f"Processed {counter} variants")
        (chrom, pos, ref, alt) = (record.CHROM, record.POS, record.REF, record.ALT[0].value)
        key = ':'.join([chrom, pos, ref, alt])
        if redis_db.exists(key):
            puts(f"variant: '{key}' already seen -- skipping")
        else:
            pipe = r.pipeline()
            variant_id = pipe.get('variant_id')
            pipe.set(key, variant_id)
            pipe.incr(key)
            pipe.sadd(redis_set, key)
            pipe.execute()
            puts(f"variant: '{key}' adding to redis -- {variant_id}")
        counter += 1

def check_global_variant_counter(redis_db):
    if not redis_db.exists('variant_id'):
        puts("Creating the variant id global counter: variant_id=1")
        redis_db.set('variant_id', 1)

def import_vcf(input_vcf, redis_host, redis_port, batch_number):
    redis_db = _redis_connect(redis_host, redis_port)
    check_global_variant_counter(redis_db)
    _process_vcf(input_vcf, redis_db, batch_number)
