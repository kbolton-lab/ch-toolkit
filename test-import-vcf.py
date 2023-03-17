#!/usr/bin/env python

# pip install
#  duckdb==0.7.0
#  ipython==8.10.0
#  vcfpy==0.13.6

import os
import vcfpy
import duckdb

vcf_file = '/scratch1/fs1/bolton/redis/playground/mutect.H_VL-MI-00052-AB43684226.vcf.gz'
variant_db_name = 'variant.db'
mutect_db_name = 'mutect.db'

def create_variants_db(db_name):
	con = duckdb.connect(db_name)

	sql = '''
    CREATE TABLE variants(
			variant_id             BIGINT PRIMARY KEY,
			chrom                  VARCHAR(5),
			pos                    INTEGER,
			ref                    varchar(255) NOT NULL,
			alt                    varchar(255) NOT NULL,
			qc_pass                boolean,
			UNIQUE( chrom, pos, ref, alt )
    )
	'''
	con.sql(sql)
	con.close()

def create_mutect_db(db_name):
	con = duckdb.connect(db_name)

	sql = '''
		CREATE TABLE IF NOT EXISTS mutect
		(
			sample_name                         varchar(25),
			variant_id                          integer,
			mutect_filter_type                  varchar(25),
			info_dp                             integer,
			info_ecnt                           integer,
			format_af                           decimal(10,5),
			format_dp                           integer,
			format_ref_count                    integer,
			format_alt_count                    integer,
			PRIMARY KEY( sample_name, variant_id )
		)
	'''
	con.sql(sql)
	con.close()

def insert_variants(vcf_file, db_name):
	con = duckdb.connect(db_name)
	reader = vcfpy.Reader.from_path(vcf_file)
	
	variant_id = 1
	for record in reader:
		chrom = record.CHROM
		pos = record.POS
		ref = record.REF
		alt = record.ALT[0].value
		qc_pass = True if record.FILTER == 'PASS' else False
		print(f"id={variant_id} | pos={pos} | ref={ref} | alt={alt} | qc_pass={qc_pass}")
		sql = f"INSERT INTO variants (variant_id, chrom, pos, ref, alt, qc_pass) VALUES ({variant_id}, '{chrom}', {pos}, '{ref}', '{alt}', {qc_pass})"
		print(sql)
		con.sql(sql)
		variant_id += 1
		if variant_id < 10:
			continue
		else:
			return

def insert_caller(vcf_file, db_name):
	con = duckdb.connect(db_name)
	reader = vcfpy.Reader.from_path(vcf_file)
	
	variant_id = 1
	sample_name = os.path.basename(vcf_file).split('.')[1]
	for record in reader:
		info_dp = record.INFO.get('DP', 0)
		format_dp = record.calls[0].data.get('DP', 0)
		mutect_filter_type = record.FILTER[0] # TODO: figure out list case
		info_ecnt = record.INFO.get('ECNT', 0)
		format_af = record.calls[0].data.get('AF')[0]
		(format_ref_count, format_alt_count) = record.calls[0].data.get('AD')
		print(f"info_dp={info_dp} | format_dp={format_dp} | mutect_filter_type={mutect_filter_type} | info_ecnt={info_ecnt} | format_af={format_af} | format_ref_count={format_ref_count} | format_alt_count={format_alt_count}")
		sql = f"INSERT INTO mutect (variant_id, sample_name, mutect_filter_type, info_dp, info_ecnt, format_af, format_dp, format_ref_count, format_alt_count) VALUES ({variant_id}, '{sample_name}', '{mutect_filter_type}', {info_dp}, {info_ecnt}, {format_af}, {format_dp}, {format_ref_count}, {format_alt_count})"
		print(sql)
		con.sql(sql)
		variant_id += 1
		if variant_id < 10:
			continue
		else:
			return

		
def select_variants(db_name):
	con = duckdb.connect(db_name)
	con.sql("SELECT * from variants where variant_id < 2").show()

def select_mutect(db_name):
	con = duckdb.connect(db_name)
	con.sql("SELECT * from mutect where variant_id < 2").show()

def join_test(variant_db_name, mutect_db_name):
	con = duckdb.connect(mutect_db_name)
	#con.sql(f"ATTACH '{variant_db_name}' (READ_ONLY)")
	join_sql = """
		select * 
		from   mutect m
		join   variant.variants v using(variant_id, variant_id)
	"""
	con.sql(join_sql).show()

# 1. variants
# 2. pon
# 3. callers (mutect, vardict)
# 4. filter types (mutect and vardict)
# 5. annotate 
# 6. sample (includes phenotypes)
# 7. fp filter
