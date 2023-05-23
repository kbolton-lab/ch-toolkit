#!/usr/bin/env python

# pip install
#  duckdb==0.7.0
#  ipython==8.10.0
#  vcfpy==0.13.6

import os
import vcfpy
import duckdb

def create_variants_tbl(db_name):
	con = duckdb.connect(db_name)

	sql = '''
	CREATE TABLE IF NOT EXISTS variants(
			variant_id             BIGINT PRIMARY KEY,
			chrom                  VARCHAR(5),
			pos                    INTEGER,
			ref                    varchar(255) NOT NULL,
			alt                    varchar(255) NOT NULL,
			qc_pass                boolean,
			batch_id               INTEGER,
			UNIQUE( chrom, pos, ref, alt )
	)
	'''
	con.execute(sql)
	con.close()

def create_samples_tbl(db_name):
	con = duckdb.connect(db_name)

	sql = '''
	CREATE TABLE IF NOT EXISTS samples(
		sample_id                  integer NOT NULL PRIMARY KEY,
		sample_name                varchar(255) NOT NULL UNIQUE
	)
	'''
	con.execute(sql)
	con.close()

def create_mutect_tbl(db_name):
	con = duckdb.connect(db_name)

	sql = '''
		CREATE TABLE IF NOT EXISTS mutect(
			sample_id                           integer NOT NULL,
			variant_id                          integer NOT NULL,
			version                             varchar(50),
			mutect_filter_type_id               varchar[],
			info_as_filterstatus                varchar[],       /* e.g. AS_FilterStatus=strand_bias|strand_bias */
			info_as_sb_table                    varchar(20),     /* e.g. AS_SB_TABLE=5379,6627|1,27|1,174 */
			info_dp                             integer,
			info_ecnt                           integer,
			info_mbq_ref                        integer,
			info_mbq_alt                        integer,
			info_mfrl_ref                       integer,
			info_mfrl_alt                       integer,
			info_mmq_ref                        integer,
			info_mmq_alt                        integer,
			info_mpos                           integer,
			info_popaf                          decimal(10,1),
			info_roq                            decimal(10,1),
			info_rpa_ref                        integer,
			info_rpa_alt                        integer,
			info_ru                             varchar(10),
			info_str                            boolean,        /* 0 Flag */
			info_strq                           integer,
			info_tlod                           decimal(10,5),
			format_af                           decimal(10,5),
			format_dp                           integer,
			format_ref_count                    integer,
			format_alt_count                    integer,
			format_ref_f1r2                     integer,
			format_alt_f1r2                     integer,
			format_ref_f2r1                     integer,
			format_alt_f2r1                     integer,
			format_gt                           varchar(3),
			format_ref_fwd                      integer,
			format_ref_rev                      integer,
			format_alt_fwd                      integer,
			format_alt_rev                      integer,
			fisher_p_value                      decimal(10,5),
			PRIMARY KEY( sample_id, variant_id )
		)
	'''
	con.execute(sql)
	con.close()

def create_mutect_filter_types_tbl(db_name):
	con = duckdb.connect(db_name)

	sql = '''
		CREATE TABLE IF NOT EXISTS mutect_filter_types(
			mutect_filter_type_id          integer NOT NULL PRIMARY KEY,
			tag                            varchar(10) NOT NULL UNIQUE
		)
	'''
	con.execute(sql)

	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (1, 'PASS');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (2, 'FAIL');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (3, 'base_qual');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (4, 'clustered_events');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (5, 'contamination');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (6, 'duplicate');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (7, 'fragment');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (8, 'germline');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (9, 'haplotype');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (10, 'low_allele_frac');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (11, 'map_qual');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (12, 'multiallelic');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (13, 'n_ratio');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (14, 'normal_artifact');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (15, 'orientation');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (16, 'panel_of_normals');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (17, 'position');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (18, 'possible_numt');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (19, 'slippage');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (20, 'strand_bias');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (21, 'strict_strand');")
	con.execute(f"INSERT INTO mutect_filter_types (mutect_filter_type_id, tag) VALUES (22, 'weak_evidence');")

	con.close()

def create_vardict_tbl(db_name):
	con = duckdb.connect(db_name)

	sql = '''
		CREATE TABLE IF NOT EXISTS vardict(
			sample_id               integer NOT NULL,
			variant_id              integer NOT NULL,
			version                 varchar(10),
			vardict_filter_type_id  varchar[],
			info_type               varchar(5),
			info_dp                 integer,
			info_vd                 integer,
			info_af                 decimal(10,5),
			info_bias               varchar(20),      /* 11671:3913 */
			info_refbias            varchar(20),
			info_varbias            varchar(20),
			info_pmean              decimal(10,5),
			info_pstd               decimal(10,5),
			info_qual               decimal(10,5),
			info_qstd               decimal(10,5),
			info_sbf                decimal(10,5),
			info_oddratio           decimal(10,5),
			info_mq                 decimal(10,5),
			info_sn                 decimal(10,5),
			info_hiaf               decimal(10,5),
			info_adjaf              decimal(10,5),
			info_shift3             integer,
			info_msi                decimal(10,5),
			info_msilen             decimal(10,5),
			info_nm                 decimal(10,5),
			info_lseq               varchar(20),
			info_rseq               varchar(20),
			info_hicnt              integer,
			info_hicov              integer,
			info_splitread          integer,
			info_spanpair           integer,
			info_duprate            decimal(10,5),
			format_ref_count        integer,
			format_alt_count        integer,
			format_gt               varchar(3),
			format_dp               integer,
			format_vd               integer,
			format_af               decimal(10,5),
			format_ref_fwd          integer,
			format_ref_rev          integer,
			format_alt_fwd          integer,
			format_alt_rev          integer,
			fisher_p_value          decimal(10,5),
			PRIMARY KEY( sample_id, variant_id )
		)
	'''
	con.execute(sql)
	con.close()

def create_vardict_filter_types_tbl(db_name):
	con = duckdb.connect(db_name)

	sql = '''
		CREATE TABLE IF NOT EXISTS vardict_filter_types(
			vardict_filter_type_id         integer NOT NULL PRIMARY KEY,
			tag                            varchar(10) NOT NULL UNIQUE
		)
	'''
	con.execute(sql)

	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (1, 'PASS');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (2, 'q22_5');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (3, 'Q10');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (4, 'p8');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (5, 'SN1_5');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (6, 'Bias');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (7, 'pSTD');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (8, 'd3');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (9, 'v2');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (10, 'min_af');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (11, 'MSI12');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (12, 'NM5_25');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (13, 'InGap');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (14, 'InIns');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (15, 'Cluster0bp');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (16, 'LongMSI');")
	con.execute(f"INSERT INTO vardict_filter_types (vardict_filter_type_id, tag) VALUES (17, 'AMPBIAS');")

	con.close()

def insert_samples(sample_name, db_name):
	con = duckdb.connect(db_name)
	sample_id = con.execute(f"SELECT MAX(sample_id) FROM samples;").fetchone()[0]
	sample_id = 1 if sample_id is None else sample_id + 1
	con.execute(f"INSERT INTO samples (sample_id, sample_name) VALUES ({sample_id}, '{sample_name}')")
	con.close()

def insert_variants(vcf_file, db_name):
	con = duckdb.connect(db_name)
	reader = vcfpy.Reader.from_path(vcf_file)

	# Check Variant ID
	variant_id = con.execute(f"SELECT MAX(variant_id) FROM variants;").fetchone()[0]
	variant_id = 0 if variant_id is None else variant_id
	batch_id = 1
	for record in reader:
		variant_id += 1
		chrom = record.CHROM
		pos = record.POS
		ref = record.REF
		alt = record.ALT[0].value
		qc_pass = True if record.FILTER == 'PASS' else False
		print(f"id={variant_id} | chrom={chrom} | pos={pos} | ref={ref} | alt={alt} | qc_pass={qc_pass} | batch_id={batch_id}")
		sql = f"INSERT INTO variants (variant_id, chrom, pos, ref, alt, qc_pass, batch_id) VALUES ({variant_id}, '{chrom}', {pos}, '{ref}', '{alt}', {qc_pass}, {batch_id})"
		print(sql)
		con.execute(sql)
		# if variant_id < 100:
		# 	continue
		# else:
		# 	return

def insert_mutect_caller(vcf_file, db_name):
	con = duckdb.connect(db_name)
	reader = vcfpy.Reader.from_path(vcf_file)
	#version = reader.header.HeaderLine("MutectVersion")    # No Idea how to get this information yet...
	version="2.2"

	sample_name = os.path.basename(vcf_file).split('.')[1]
	for record in reader:
		# Get the Variant_ID
		variant_id = con.execute(f"SELECT variant_id FROM variants WHERE chrom = '{record.CHROM}' AND pos = {record.POS} AND ref = '{record.REF}' AND alt = '{record.ALT[0].value}';").fetchone()[0]
		# Get the Sample_ID
		sample_id = con.execute(f"SELECT sample_id FROM samples WHERE sample_name = '{sample_name}';").fetchone()[0]

		mutect_filter_type = record.FILTER
		info_as_filterstatus = ['Multiallelic'] if record.INFO.get('AS_FilterStatus') is None else record.INFO.get('AS_FilterStatus')
		info_as_sb_table = record.INFO.get('AS_SB_TABLE')
		info_dp = record.INFO.get('DP', 0)
		info_ecnt = record.INFO.get('ECNT', 0)
		info_mbq_ref = record.INFO.get("MBQ")[0]
		info_mbq_alt = record.INFO.get("MBQ")[1]
		info_mfrl_ref = record.INFO.get("MFRL")[0]
		info_mfrl_alt = record.INFO.get("MFRL")[1]
		info_mmq_ref = record.INFO.get("MMQ")[0]
		info_mmq_alt = record.INFO.get("MMQ")[1]
		info_mpos = record.INFO.get("MPOS")[0]
		info_popaf = record.INFO.get("POPAF")[0]
		info_roq = record.INFO.get("ROQ")
		info_rpa_ref = 'NULL' if record.INFO.get("RPA") is None else record.INFO.get("RPA")[0]
		info_rpa_alt = 'NULL' if record.INFO.get("RPA") is None else record.INFO.get("RPA")[1]
		info_ru = 'NULL' if record.INFO.get("RU") is None else record.INFO.get("RU")[0]
		info_str = record.INFO.get("STR", False)
		info_strq = record.INFO.get("STRQ", "NULL")
		info_tlod = record.INFO.get("TLOD")[0]

		format_dp = record.calls[0].data.get('DP', 0)
		format_af = record.calls[0].data.get('AF')[0]
		(format_ref_count, format_alt_count) = record.calls[0].data.get('AD')
		(format_ref_f1r2, format_alt_f1r2) = record.calls[0].data.get('F1R2')
		(format_ref_f2r1, format_alt_f2r1) = record.calls[0].data.get('F2R1')
		format_gt = record.calls[0].data.get('GT')
		(format_ref_fwd, format_ref_rev, format_alt_fwd, format_alt_rev) = record.calls[0].data.get('SB')
		fisher_p_value = 0.05

		print(f"sample_id={sample_id} | variant_id={variant_id} | version={version} | mutect_filter_type={mutect_filter_type} \nINFO: info_as_filterstatus={info_as_filterstatus} | info_as_sb_table={info_as_sb_table} | info_dp={info_dp} | info_ecnt={info_ecnt} | info_mbq_ref={info_mbq_ref} | info_mbq_alt={info_mbq_alt} | info_mfrl_ref={info_mfrl_ref} | info_mfrl_alt={info_mfrl_alt} | info_mmq_ref={info_mmq_ref} | info_mmq_alt={info_mmq_alt} | info_mpos={info_mpos} | info_popaf={info_popaf} | info_roq={info_roq} | info_rpa_ref={info_rpa_ref} | info_rpa_alt={info_rpa_alt} | info_ru={info_ru} | info_str={info_str} | info_strq={info_strq} | info_tlod={info_tlod} \nFORMAT: format_dp={format_dp} | format_af={format_af} | format_ref_count={format_ref_count} | format_alt_count={format_alt_count} | format_ref_f1r2={format_ref_f1r2} | format_alt_f1r2={format_alt_f1r2} | format_ref_f2r1={format_ref_f2r1} | format_alt_f2r1={format_alt_f2r1} | format_gt={format_gt} | format_ref_fwd={format_ref_fwd} | format_ref_rev={format_ref_rev} | format_alt_fwd={format_alt_fwd} | format_alt_rev={format_alt_rev} | fisher_p_value={fisher_p_value}")
		sql = f"INSERT INTO mutect (sample_id, variant_id, version, mutect_filter_type_id, info_as_filterstatus, info_as_sb_table, info_dp, info_ecnt, info_mbq_ref, info_mbq_alt, info_mfrl_ref, info_mfrl_alt, info_mmq_ref, info_mmq_alt, info_mpos, info_popaf, info_roq, info_rpa_ref, info_rpa_alt, info_ru, info_str, info_strq, info_tlod, format_gt, format_af, format_dp, format_ref_count, format_alt_count, format_ref_f1r2, format_alt_f1r2, format_ref_f2r1, format_alt_f2r1, format_ref_fwd, format_ref_rev, format_alt_fwd, format_alt_rev, fisher_p_value) VALUES ({sample_id}, {variant_id}, '{version}', {mutect_filter_type}, {info_as_filterstatus}, '{info_as_sb_table}', {info_dp}, {info_ecnt}, {info_mbq_ref}, {info_mbq_alt}, {info_mfrl_ref}, {info_mfrl_alt}, {info_mmq_ref}, {info_mmq_alt}, {info_mpos}, {info_popaf}, {info_roq}, {info_rpa_ref}, {info_rpa_alt}, '{info_ru}', {info_str}, {info_strq}, {info_tlod}, {format_gt}, {format_af}, {format_dp}, {format_ref_count}, {format_alt_count}, {format_ref_f1r2}, {format_alt_f1r2}, {format_ref_f2r1}, {format_alt_f2r1}, {format_ref_fwd}, {format_ref_rev}, {format_alt_fwd}, {format_alt_rev}, {fisher_p_value})"
		print(sql + "\n")
		con.execute(sql)

		# if variant_id < 100:
		# 	continue
		# else:
		# 	return

def insert_vardict_caller(vcf_file, db_name):
	con = duckdb.connect(db_name)
	reader = vcfpy.Reader.from_path(vcf_file)
	version="v1.8.2"

	sample_name = os.path.basename(vcf_file).split('.')[1]
	for record in reader:
		# Get the Variant_ID
		variant_id = con.execute(f"SELECT variant_id FROM variants WHERE chrom = '{record.CHROM}' AND pos = {record.POS} AND ref = '{record.REF}' AND alt = '{record.ALT[0].value}';").fetchone()[0]
		# Get the Sample_ID
		sample_id = con.execute(f"SELECT sample_id FROM samples WHERE sample_name = '{sample_name}';").fetchone()[0]

		filter_type = record.FILTER
		info_type = record.INFO.get('TYPE')
		info_dp = record.INFO.get('DP', 0)
		info_vd = record.INFO.get('VD', 0)
		info_af = record.INFO.get('AF', 0)[0]
		info_bias = record.INFO.get('BIAS')
		info_refbias = record.INFO.get('REFBIAS')
		info_varbias = record.INFO.get('VARBIAS')
		info_pmean = record.INFO.get('PMEAN')
		info_pstd = record.INFO.get('PSTD')
		info_qual = record.INFO.get('QUAL')
		info_qstd = record.INFO.get('QSTD')
		info_sbf = record.INFO.get('SBF')
		info_oddratio = record.INFO.get('ODDRATIO')
		info_mq = record.INFO.get('MQ')
		info_sn = record.INFO.get('SN')
		info_hiaf = record.INFO.get('HIAF')
		info_adjaf = record.INFO.get('ADJAF')
		info_shift3 = record.INFO.get('SHIFT3')
		info_msi = record.INFO.get('MSI')
		info_msilen = record.INFO.get('MSILEN')
		info_nm = record.INFO.get('NM')
		info_lseq = record.INFO.get('LSEQ')
		info_rseq = record.INFO.get('RSEQ')
		info_hicnt = record.INFO.get('HICNT')
		info_hicov = record.INFO.get('HICOV')
		info_splitread = record.INFO.get('SPLITREAD')
		info_spanpair = record.INFO.get('SPANPAIR')
		info_duprate = record.INFO.get('DUPRATE')

		format_gt = record.calls[0].data.get('GT')
		format_dp = record.calls[0].data.get('DP', 0)
		format_af = record.calls[0].data.get('AF')[0]
		(format_ref_count, format_alt_count) = record.calls[0].data.get('AD')
		(format_ref_fwd, format_ref_rev) = record.calls[0].data.get('RD')
		(format_alt_fwd, format_alt_rev) = record.calls[0].data.get('ALD')
		format_vd = record.calls[0].data.get('VD', 0)
		fisher_p_value = 0.05

		print(f"sample_id={sample_id} | variant_id={variant_id} | version={version} | vardict_filter_type={filter_type} \nINFO: info_type={info_type} | info_vd={info_vd} | info_af={info_af} | info_bias={info_bias} | info_refbias={info_refbias} | info_varbias={info_varbias} | info_pmean={info_pmean} | info_pstd={info_pstd} | info_qual={info_qual} | info_qstd={info_qstd} | info_sbf={info_sbf} | info_oddratio={info_oddratio} | info_mq={info_mq} | info_sn={info_sn} | info_hiaf={info_hiaf} | info_adjaf={info_adjaf} | info_shift3={info_shift3} | info_msi={info_msi} | info_msilen={info_msilen} | info_nm={info_nm} | info_lseq={info_lseq} | info_rseq={info_rseq} | info_hicnt={info_hicnt} | info_hicov={info_hicov} | info_splitread={info_splitread} | info_spanpair={info_spanpair} | info_duprate={info_duprate} \nFORMAT: format_dp={format_dp} | format_af={format_af} | format_ref_count={format_ref_count} | format_alt_count={format_alt_count} | format_vd={format_vd} | format_gt={format_gt} | format_ref_fwd={format_ref_fwd} | format_ref_rev={format_ref_rev} | format_alt_fwd={format_alt_fwd} | format_alt_rev={format_alt_rev} | fisher_p_value={fisher_p_value}")
		sql = f"INSERT INTO vardict (sample_id, variant_id, version, vardict_filter_type_id, info_type, info_vd, info_af, info_bias, info_refbias, info_varbias, info_pmean, info_pstd, info_qual, info_qstd, info_sbf, info_oddratio, info_mq, info_sn, info_hiaf, info_adjaf, info_shift3, info_msi, info_msilen, info_nm, info_lseq, info_rseq, info_hicnt, info_hicov, info_splitread, info_spanpair, info_duprate, format_gt, format_af, format_dp, format_ref_count, format_alt_count, format_vd, format_ref_fwd, format_ref_rev, format_alt_fwd, format_alt_rev, fisher_p_value) VALUES ({sample_id}, {variant_id}, '{version}', {filter_type}, '{info_type}', {info_vd}, {info_af}, '{info_bias}', '{info_refbias}', '{info_varbias}', {info_pmean}, {info_pstd}, {info_qual}, {info_qstd}, {info_sbf}, {info_oddratio}, {info_mq}, {info_sn}, {info_hiaf}, {info_adjaf}, {info_shift3}, {info_msi}, {info_msilen}, {info_nm}, '{info_lseq}', '{info_rseq}', {info_hicnt}, {info_hicov}, {info_splitread}, {info_spanpair}, {info_duprate}, {format_gt}, {format_af}, {format_dp}, {format_ref_count}, {format_alt_count}, {format_vd}, {format_ref_fwd}, {format_ref_rev}, {format_alt_fwd}, {format_alt_rev}, {fisher_p_value})"
		print(sql + "\n")
		con.execute(sql)

		# if variant_id < 100:
		# 	continue
		# else:
		# 	return

def insert_caller(vcf_file, db_name):
	caller_name = os.path.basename(vcf_file).split('.')[0]
	if caller_name == "mutect":
		insert_mutect_caller(vcf_file, db_name)
	elif caller_name == "vardict":
		insert_vardict_caller(vcf_file, db_name)

def select_variants(db_name):
	con = duckdb.connect(db_name)
	con.execute("SELECT * from variants where variant_id < 2").show()

def select_mutect(db_name):
	con = duckdb.connect(db_name)
	con.execute("SELECT * from mutect where variant_id < 2").show()

def join_test(variant_db_name, mutect_db_name):
	con = duckdb.connect(mutect_db_name)
	join_sql = """
		select *
		from   mutect m
		join   variant.variants v using(variant_id, variant_id)
	"""
	con.execute(join_sql).show()

if __name__ == '__main__':
	vcf_file_mutect = '/storage1/fs1/bolton/Active/Projects/chip-toolkit/data/TERRA_Data/H_VL-MI-00052-AB43684226/mutect.H_VL-MI-00052-AB43684226.vcf.gz'
	vcf_file_vardict = '/storage1/fs1/bolton/Active/Projects/chip-toolkit/data/TERRA_Data/H_VL-MI-00052-AB43684226/vardict.H_VL-MI-00052-AB43684226.vcf.gz'
	variant_db_name = 'variant.db'
	mutect_db_name = 'mutect.db'
	project_db = 'ccdg_v.db'

	con = duckdb.connect(project_db)

	create_samples_tbl(project_db)
	create_variants_tbl(project_db)
	create_mutect_tbl(project_db)
	create_mutect_filter_types_tbl(project_db)
	create_vardict_tbl(project_db)
	create_vardict_filter_types_tbl(project_db)

	insert_samples("H_VL-MI-00052-AB43684226", project_db)

	#insert_variants(vcf_file_mutect, project_db)
	insert_variants(vcf_file_vardict, project_db)

	#insert_caller(vcf_file_mutect, project_db)
	insert_caller(vcf_file_vardict, project_db)

	con.close()


# 1. variants
# 2. pon
# 3. callers (mutect, vardict)
# 4. filter types (mutect and vardict)
# 5. annotate
# 6. sample (includes phenotypes)
# 7. fp filter
