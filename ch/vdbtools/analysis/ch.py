import sys
import math
import duckdb
import pandas as pd
import ch.utils.logger as log
import ch.utils.database as db
import importlib.resources
import numpy as np

def load_flat_databases():
    bolton_bick_vars = importlib.resources.files('ch.resources.annotate_pd').joinpath('bick.bolton.vars3.txt')
    cosmic_hotspots = importlib.resources.files('ch.resources.annotate_pd').joinpath('COSMIC.heme.myeloid.hotspot.w_truncating_counts.tsv')
    pd_table = importlib.resources.files('ch.resources.annotate_pd').joinpath('pd_table_kbreview_bick_trunc4_oncoKB_SAFE.filtered_genes_oncoKB_CGC.tsv')
    gene_list = importlib.resources.files('ch.resources.annotate_pd').joinpath('oncoKB_CGC_pd_table_disparity_KB_BW.csv')

    vars = pd.read_csv(bolton_bick_vars, sep='\t')
    vars['aa.pos'] = vars['loci.vep'].str.extract(r'(\d+)').astype(float)                               # e.g. DNMT3A_R882 --> 882
    vars['CHROM.POS'] = vars['CHROM'] + '_' + vars['POS'].astype(str)                                   # e.g. chr2_25234373
    vars['GENE.AA.POS'] = vars.apply(lambda x: x['SYMBOL_VEP'] + '_' + str(int(x['aa.pos'])) if not pd.isna(x['aa.pos']) else pd.NA, axis=1)    # e.g. DNMT3A_882
    vars.loc[vars['GENE.AA.POS'] == 'NA_NA', 'GENE.AA.POS'] = pd.NA
    vars['gene_cDNAchange'] = vars['SYMBOL_VEP'] + '_' + vars['HGVSc_VEP'].str.extract(r'(c\.\d+[A-Z]>[A-Z])', expand=False).iloc[0]       # e.g. DNMT3A_c.2645G>A
    vars['gene_aachange'] = vars['SYMBOL_VEP'] + '_' + vars['AAchange2']                                # e.g. DNMT3A_R882H

    ct = pd.read_csv(cosmic_hotspots, sep='\t')
    ct['CHROM.POS'] = ct['CHROM'] + '_' + ct['POS'].astype(str)
    ct['GENE.AA.POS'] = ct['gene'] + '_' + ct['aa.pos'].astype(str)
    ct['gene_cDNAchange'] = ct['gene'] + '_' + ct['cDNAchange']
    ct['gene_aachange'] = ct['gene'] + '_' + vars['AAchange']

    bickGene = pd.read_csv(pd_table, sep='\t')
    bickGene = bickGene[bickGene['source'] == 'Bick_email']
    gene_list = pd.read_csv(gene_list, sep=',')
    TSG_gene_list = gene_list[gene_list['isTSG'] == 1]['Gene']
    gene_list = gene_list['Gene']

    ZBTB33 = bickGene.loc[bickGene['Gene'] == "ZBTB33", ['aa_ref', 'aa_pos', 'aa_alt']]
    ZBTB33 = ZBTB33.loc[ZBTB33['aa_ref'] != "***", ['aa_ref', 'aa_pos', 'aa_alt']]
    ZBTB33['AAchange'] = ZBTB33['aa_ref'] + ZBTB33['aa_pos'] + ZBTB33['aa_alt']
    ZBTB33 = ZBTB33['AAchange'].unique()

    return vars, ct, bickGene, TSG_gene_list, gene_list, ZBTB33
    
def near_BB_loci_HS(df_row, vars):
    p = list(range(-3, 4))              # Establishes 3 AA downstream and upstream
    n = list(range(-9, 10))             # Establishes the 9 nucleotides downstream and upstream

    prot = [str(int(df_row['aa.pos'] + x)) for x in p]                  # e.g. DNMT3A_R882 would result in [879, 880, 881, 882, 883, 884, 885]
    vector_p = pd.Series([df_row['SYMBOL'] + '_' + x for x in prot])    # which would result in [DNMT3A_879, DNMT3A_880, ...]
    any_in_p = vector_p.isin(vars['GENE.AA.POS'])                       # If any of these AA positions are inside vars, it means this mutation is NEAR a BB Hotspot
    
    nuc = [str(int(df_row['POS']) + x) for x in n]                      # e.g. DNMT3A_R882 would result in [25234364, 25234365, ...]
    vector_n = pd.Series([df_row['CHROM'] + '_' + x for x in nuc])      # which would result in [chr2_25234364, chr2_25234365, ...]
    any_in_n = vector_n.isin(vars['CHROM.POS'])                         # Again, if any of these specific positions are inside vars... it means it is NEAR a BB Hotspot

    # If at least one of the AA POS is true...
    if any(any_in_p):  
        # Grabs the rows from vars as long as it is a Hotspot e.g. DNMT3A_R882 0 663
        res = vars[(vars['GENE.AA.POS'].isin(vector_p)) & (vars['n.loci.truncating.vep'] >= 5)][['gene_loci_vep', 'truncating', 'n.loci.truncating.vep']].drop_duplicates()        
        # If the AA change is Terminating, we only want to select the Near Hotspots that are ALSO Terminating, or else it's comparing Apples to Oranges
        if 'Ter' in df_row['gene_aachange']:
            res = res[res['truncating'] == 'truncating']
            return res.apply(lambda row: f"{row['gene_loci_vep']}: {row['n.loci.truncating.vep'].astype(str)}", axis=1).str.cat(sep=' | ')
            #return res.apply(lambda row: f"{row['gene_loci_vep']}: {row['n.loci.truncating.vep']}", axis=1).str.cat(sep=' | ')
        else:
            res = res[res['truncating'] == 'not']
            return res.apply(lambda row: f"{row['gene_loci_vep']}: {row['n.loci.truncating.vep'].astype(str)}", axis=1).str.cat(sep=' | ')
            #return res.apply(lambda row: f"{row['gene_loci_vep']}: {row['n.loci.truncating.vep']}", axis=1).str.cat(sep=' | ')
    # Else check nucleotide if none of the AA positions are true...
    elif any(any_in_n):
        res = vars[vars['CHROM.POS'].isin(vector_n)][['gene_loci_vep', 'n.HGVSc', 'gene_cDNAchange']].drop_duplicates()
        # Similar to Terminating, if the cDNAchange is a deletion, insertion, or duplication, we only want to compare to those results
        if any(res['gene_cDNAchange'].str.contains('del|ins|dup')):
            res = res[res['gene_cDNAchange'].str.contains('del|ins|dup')]
            res = res.groupby(['n.HGVSc', 'gene_loci_vep']).size().reset_index(name='n')
            return res.apply(lambda row: f"{row['gene_loci_vep']}: {row['n'].astype(str)}", axis=1).str.cat(sep=' | ')
        else:
            res = res[~res['gene_cDNAchange'].str.contains('del|ins|dup')]
            res = res.groupby('gene_loci_vep').size().reset_index(name='n')
            return res.apply(lambda row: f"{row['gene_loci_vep']}: {row['n'].astype(str)}", axis=1).str.cat(sep=' | ')
    else:
        return ''
    
def near_COSMIC_loci_HS(df_row, ct):
    p = list(range(-3, 4))              # Establishes 3 AA downstream and upstream
    n = list(range(-9, 10))             # Establishes the 9 nucleotides downstream and upstream
    
    prot = [str(int(df_row['aa.pos'] + x)) for x in p]                  # e.g. DNMT3A_R882 would result in [879, 880, 881, 882, 883, 884, 885]
    vector_p = pd.Series([df_row['SYMBOL'] + '_' + x for x in prot])    # which would result in [DNMT3A_879, DNMT3A_880, ...]
    any_in_p = vector_p.isin(ct['GENE.AA.POS'])                       # If any of these AA positions are inside COSMIC, it means this mutation is NEAR a BB Hotspot
    
    nuc = [str(int(df_row['POS']) + x) for x in n]                      # e.g. DNMT3A_R882 would result in [25234364, 25234365, ...]
    vector_n = pd.Series([df_row['CHROM'] + '_' + x for x in nuc])      # which would result in [chr2_25234364, chr2_25234365, ...]
    any_in_n = vector_n.isin(ct['CHROM.POS'])                         # Again, if any of these specific positions are inside COSMIC... it means it is NEAR a BB Hotspot
    
    # If at least one of the AA POS is true...
    if any(any_in_p):  
        res = ct[(ct['GENE.AA.POS'].isin(vector_p)) & ((ct['cosmic_count.loci.truncating'] >= 25) | (ct['heme_count.loci.truncating'] >= 10) | (ct['myeloid_count.loci.truncating'] >= 5))][['gene_loci_vep', 'truncating', 'cosmic_count.loci.truncating', 'heme_count.loci.truncating', 'myeloid_count.loci.truncating']].drop_duplicates()
        # If the AA change is Terminating, we only want to select the Near Hotspots that are ALSO Terminating, or else it's comparing Apples to Oranges
        if 'Ter' in df_row['gene_aachange']:
            res = res[res['truncating'] == True]
            return res.apply(lambda row: f"{row['gene_loci_vep']}: {row['cosmic_count.loci.truncating']}, {row['heme_count.loci.truncating']}, {row['myeloid_count.loci.truncating']}", axis=1).str.cat(sep=' | ')
        else:
            res = res[res['truncating'] == False]
            return res.apply(lambda row: f"{row['gene_loci_vep']}: {row['cosmic_count.loci.truncating']}, {row['heme_count.loci.truncating']}, {row['myeloid_count.loci.truncating']}", axis=1).str.cat(sep=' | ')
    # Else check nucleotide if none of the AA positions are true...
    elif any(any_in_n):
        res = ct[ct['CHROM.POS'].isin(vector_n)][['gene_loci_vep', 'cosmic_count.totals.c', 'heme_count.totals.c', 'myeloid_count.totals.c', 'gene_cDNAchange']].drop_duplicates()
        # Similar to Terminating, if the cDNAchange is a deletion, insertion, or duplication, we only want to compare to those results
        if any(res['gene_cDNAchange'].str.contains('del|ins|dup')):
            res = res[res['gene_cDNAchange'].str.contains('del|ins|dup')]
            res = res.groupby('gene_loci_vep')[['cosmic_count.totals.c', 'heme_count.totals.c', 'myeloid_count.totals.c']].sum().reset_index()
            res = res[(res['cosmic_count.totals.c'] >= 25) | (res['heme_count.totals.c'] >= 10) | (res['myeloid_count.totals.c'] >= 5)]
            return res.apply(lambda row: f"{row['gene_loci_vep']}: {row['cosmic_count.totals.c']}, {row['heme_count.totals.c']}, {row['myeloid_count.totals.c']}", axis=1).str.cat(sep=' | ')
        else:
            res = res[~res['gene_cDNAchange'].str.contains('del|ins|dup')]
            res = res.groupby('gene_loci_vep')[['cosmic_count.totals.c', 'heme_count.totals.c', 'myeloid_count.totals.c']].sum().reset_index()
            res = res[(res['cosmic_count.totals.c'] >= 25) | (res['heme_count.totals.c'] >= 10) | (res['myeloid_count.totals.c'] >= 5)]
            return res.apply(lambda row: f"{row['gene_loci_vep']}: {row['cosmic_count.totals.c']}, {row['heme_count.totals.c']}, {row['myeloid_count.totals.c']}", axis=1).str.cat(sep=' | ')
    else:
        return ''

def determine_pathogenicity(df, total_samples, debug):
    # Autofail Variants with NO support in B/B and Cosmic that have nsamples >20 and in non complex region
    vars, ct, bickGene, TSG_gene_list, gene_list, ZBTB33 = load_flat_databases()

    df['sample_key'] = df['sample_name'] + df['key']
    df = df[~(df['n.HGVSc'].isna() & (df['heme_cosmic_count'] == 0) & (df['n_samples'] > 20) & (df['homopolymerCase'] == "") & (df['dust_score'] <= 7) & (~df['Gene'].isin(bickGene['Gene'])))]
    
    # Autofail Variants that are too recurrent
    bbCutoff = max(math.ceil(total_samples * 0.005), 3) # Round of total samples * 0.5% or 3
    if debug: log.logit(f"n_samples cutoff is: {bbCutoff}")
    df['n.HGVSc'] = df['n.HGVSc'].fillna(0)
    df['n.HGVSp'] = df['n.HGVSp'].fillna(0)
    df['n.loci.vep'] = df['n.loci.vep'].fillna(0)
    df['n.loci.truncating.vep'] = df['n.loci.truncating.vep'].fillna(0)
    df['n.HGVSc'] = pd.to_numeric(df['n.HGVSc'], errors='coerce')
    df['n.HGVSp'] = pd.to_numeric(df['n.HGVSp'], errors='coerce')
    df['n.loci.vep'] = pd.to_numeric(df['n.loci.vep'], errors='coerce')
    df['n.loci.truncating.vep'] = pd.to_numeric(df['n.loci.truncating.vep'], errors='coerce')
    df = df[(df['n_samples'] < bbCutoff) | (df['n.HGVSc'] > 25) | (df['CosmicCount'] > 50)]

    df[['CHROM', 'POS', 'REF', 'ALT']] = df['key'].str.split(':', expand=True)
    new_order = ['CHROM', 'POS', 'REF', 'ALT'] + [col for col in df.columns if col not in ['CHROM', 'POS', 'REF', 'ALT']]
    df = df.reindex(columns=new_order)
    df = df[(df['REF'].str.len() <= 100) & (df['ALT'].str.len() <= 100)]                                    # long100_indel
    df = df[((df['REF'].str.len() <= 20) & (df['ALT'].str.len() <= 20)) | (df['CosmicCount'] > 0)]          # long_indel
    df = df[(df['REF'].str.len() != df['ALT'].str.len()) | (df['ALT'].str.len() <= 1)]                      # di_tri_nuc

    # Adding Near Hotspots
    df['aa.pos'] = df['AAchange'].str.extract(r'(\d+)').astype(float)
    df['near.BB.loci.HS'] = df.apply(lambda row: near_BB_loci_HS(row, vars) if not pd.isna(row['aa.pos']) else "", axis=1)
    df['near.COSMIC.loci.HS'] = df.apply(lambda row: near_COSMIC_loci_HS(row, ct) if not pd.isna(row['aa.pos']) else "", axis=1)
    df['nearBBLogic'] = df['near.BB.loci.HS'] == ''
    df['nearCosmicHemeLogic'] = df['near.COSMIC.loci.HS'] == ''
    
    # Recurrent Proportion
    #total_samples = 48212           # This is inside the CCDG
    total_heme_cosmic = 29234
    total_BB = 122691

    df['n.loci.truncating.vep'] = df['n.loci.truncating.vep'].fillna(0)
    df['prop_nsamples'] = df['n_samples'] / total_samples
    df['prop_Cosmic'] = (0.5 + df['heme_cosmic_count']) / total_heme_cosmic
    df['prop_BB'] = (0.5 + df['n.HGVSc']) / total_BB
    df['ratio_to_BB'] = df['prop_nsamples'] / df['prop_BB']
    df['ratio_to_cosmic'] = df['prop_nsamples'] / df['prop_Cosmic']

    df['pass_prop_recurrent'] = np.where(
        df['Gene'].isin(bickGene['Gene']), True,
        np.where(
            (df['n_samples'] > 5) & (df['n.HGVSc'] == 0) & (df['heme_cosmic_count'] != 0) & (df['ratio_to_cosmic'] > 10), False,
            np.where(
                (df['n_samples'] > 5) & (df['n.HGVSc'] != 0) & (df['heme_cosmic_count'] == 0) & (df['ratio_to_BB'] > 6), False,
                np.where(
                    (df['n_samples'] > 5) & (df['n.HGVSc'] != 0) & (df['heme_cosmic_count'] != 0) & ((df['ratio_to_BB'] > 6) & (df['ratio_to_cosmic'] > 10)), False, True
                )
            )
        )
    )
    # df = df[(data['pass_prop_recurrent']) | ((df['key'].isin(["chr20:32434638:A:AG", "chr20:32434638:A:AGG"])) & (df['average_af'] >= 0.05))]

    # Basic Definitions
    nonsense_mutation = pd.Series(["frameshift_variant", "stop_lost", "stop_gained", "transcript_ablation"])
    missense_mutation = pd.Series(["missense_variant", "inframe_deletion", "inframe_insertion"])
    SF3B1_positions = pd.Series([622, 623, 624, 625, 626, 662, 663, 664, 665, 666, 700, 701, 702, 703, 704, 740, 741, 742])
    clinvar_sig_terms = pd.Series(["Likely_pathogenic", "Pathogenic", "Pathogenic/Likely_pathogenic", "Pathogenic/Likely_pathogenic|risk_factor", "Pathogenic|drug_response|other"])
    splicingSynonymous = pd.Series([",".join(["splice_donor_5th_base_variant", "splice_donor_region_variant", "splice_region_variant", "splice_polypyrimidine_tract_variant", "synonymous_variant"])])

    # Putative Driver Rules
    # ---------------------
    # 1) TSG + Nonsense Mutation --> PD = 1
    # 2) OncoKB is Reviewed by Pathologists --> PD = 1
    # 3) If OncoKB Reports as 'Neutral' but A LOT of Support from B/B then B/B takes precedence
    # 4) OncoKB No Support --> PD = 0
    # 5) Missense Variant + Cosmic Support  --> PD = 1
    # 6) Missense Variant + B/B Loci Count + SIFT & PolyPhen Support --> PD = 1
    # 7) Missense Variant + Near B/B Hotspot | Near Cosmic Hotspot + SIFT & PolyPhen Support --> PD = 1
    # 8a) SRSF2 + Hotspot
    # 8b) SRSF2 + OncoKB
    # 9) SF3B1 Rules
    # 10) IDH1 and IDH2 Hotspots
    # 11) JAK2 Rules
    # 12) PPM1D Exon 6 Rules
    # 13) Missense Variant + B/B AA Support w/ EITHER SIFT | PolyPhen Support --> PD = 1
    # 14) Missense Variant + Extension of Termination Codon --> PD = 1
    # 15) TSG + Splice Acceptor/Donor Variant --> PD = 1
    # 16) TSG + ClinVar Support --> PD == 1
    # 17) Synonymous Variant in Splicing Region (Handled with SpliceAI)
    # 18) Bick's Email Rules - ZBTB33
    # 19) Bick's Email Rules - Other Genes
    df['pd_reason'] = np.select(
        [
            df['Gene'].isin(TSG_gene_list) & df['VariantClass'].isin(nonsense_mutation),
            df['oncoKB'].str.contains("Oncogenic"),
            df['Gene'].isin(gene_list) & df['VariantClass'].isin(missense_mutation) & ((df['n.HGVSp'] >= 10) | (df['n.HGVSc'] >= 5)),
            df['Gene'].isin(gene_list) & df['oncoKB'].str.contains("Neutral"),
            df['Gene'].isin(gene_list) & df['VariantClass'].isin(missense_mutation) & ((df['CosmicCount'] >= 10) | (df['heme_cosmic_count'] >= 5) | (df['myeloid_cosmic_count'] >= 1)),
            df['Gene'].isin(gene_list) & df['VariantClass'].isin(missense_mutation) & ((df['n.loci.vep'] - df['n.loci.truncating.vep']) >= 5) & df['SIFT'].str.contains("deleterious") & df['PolyPhen'].str.contains("damaging"),
            df['Gene'].isin(gene_list) & df['VariantClass'].isin(missense_mutation) & ((df['nearBBLogic']) | (df['nearCosmicHemeLogic'])) & df['SIFT'].str.contains("deleterious") & df['PolyPhen'].str.contains("damaging"),
            (df['Gene'] == "SRSF2") & df['VariantClass'].isin(missense_mutation) & (df['aa.pos'] == 95),
            (df['Gene'] == "SF3B1") & df['VariantClass'].isin(missense_mutation) & df['aa.pos'].isin(SF3B1_positions),
            (df['Gene'] == "IDH1") & df['VariantClass'].isin(missense_mutation) & (df['aa.pos'] == 132),
            (df['Gene'] == "IDH2") & df['VariantClass'].isin(missense_mutation) & ((df['aa.pos'] == 140) | (df['aa.pos'] == 172)),
            (df['Gene'] == "PPM1D") & df['VariantClass'].isin(nonsense_mutation) & (df['EXON'] == "6/6"),
            df['Gene'].isin(gene_list) & df['VariantClass'].isin(missense_mutation) & (df['n.HGVSp'] >= 1) & (df['SIFT'].str.contains("deleterious") | df['PolyPhen'].str.contains("damaging")),
            df['Gene'].isin(TSG_gene_list) & df['VariantClass'].isin(["splice_donor_variant", "splice_aceptor_variant", "splice_region_variant"]) & ~(df['Consequence'].isin(splicingSynonymous)),
            df['Gene'].isin(TSG_gene_list) & df['clinvar_CLNSIG'].isin(clinvar_sig_terms),
            df['Gene'].isin(gene_list) & df['Consequence'].isin(splicingSynonymous),
            (df['Gene'] == "ZBTB33") & df['VariantClass'].isin(missense_mutation) & df['AAchange'].isin(ZBTB33),
            df['Gene'].isin(bickGene['Gene'].unique()) & df['VariantClass'].isin(nonsense_mutation)
        ],
        [
            "Nonsense Mutation in TSG",
            "OncoKB",
            "B/B Hotspot >= 10",
            "Not PD",
            "COSMIC",
            "Loci + SIFT/PolyPhen",
            "Near Hotspot + SIFT/PolyPhen",
            "SRSF2 Hotspot",
            "SF3B1 Hotspot",
            "IDH1 Hotspot",
            "IDH2 Hotspot",
            "Nonsense Mutation on PPM1D Exon 6",
            "B/B Hotspot + SIFT/PolyPhen",
            "Splicing Mutation",
            "ClinVar",
            "Not PD",
            "Bick's Email",
            "Bick's Email"
        ],
        default="Not PD"
    )
    # OncoKB API automatically classifies ALL splicing mutations as oncogenic, however if the mutation is synonymous, then we have to change it
    df.loc[(df['oncoKB'].str.contains("Oncogenic")) & (df['Consequence'].isin(splicingSynonymous)) & ~(df['clinvar_CLNSIG'].isin(clinvar_sig_terms)), 'pd_reason'] = "Not PD"
    df['putative_driver'] = np.where(df['pd_reason'] != "Not PD", 1, 0)

    # SpliceAI
    if 'SpliceAI_pred' in df.columns:
        df[['SpliceAI_pred_SYMBOL', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_DP_AG', 'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL']] = df['SpliceAI_pred'].str.split('|', expand=True)
        df[['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_DP_AG', 'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL']] = df[['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_DP_AG', 'SpliceAI_pred_DP_AL', 'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL']].apply(pd.to_numeric)

    # Start of Review
    df['Review'] = "No Review"

    # Split Protein_position into beginning and end
    df['Protein_position_start'] = pd.to_numeric(np.where(df['Protein_position'].str.contains('-'), df['Protein_position'].str.split('-').str[0], df['Protein_position']), errors='coerce')
    df['Protein_position_end'] = pd.to_numeric(np.where(df['Protein_position'].str.contains('-'), df['Protein_position'].str.split('-').str[1], df['Protein_position']), errors='coerce')

    # If it is an Inframe Insertion or Deletion that OVERLAPS P95, then mark for review
    df['Review'] = np.where((df['Review'] == "No Review") & (df['Gene'] == "SRSF2") & (df['VariantClass'].isin(["inframe_deletion", "inframe_insertion"])) & (df['Protein_position_start'].astype(float) <= 95) & (df['Protein_position_end'].astype(float) >= 95), "MW Review", df['Review'])

    # For SRSF2, SF3B1, IDH1, IDH2, and JAK2 (If not already identified as PD or review... remove all other variants)
    SSIIJ = pd.Series(["SRSF2", "SF3B1", "IDH1", "IDH2", "JAK2"])
    df = df[~((df['Gene'].isin(SSIIJ)) & (df['putative_driver'] == 0) & (df['Review'] == "No Review"))]

    unique_genes = bickGene[bickGene['Gene'] != "ZBTB33"]['Gene'].unique()
    conditions = [
        ((df['Review'] == "No Review") & (df['REF'].str.len() > 5) | (df['ALT'].str.len() > 5), "Long INDEL"),
        ((df['Review'] != "No Review") & (df['REF'].str.len() > 5) | (df['ALT'].str.len() > 5), ";Long INDEL"),
        ((df['Review'] == "No Review") & (df['REF'].str.len() >= 2) & (df['ALT'].str.len() >= 2), "Complex INDEL"),
        ((df['Review'] != "No Review") & (df['REF'].str.len() >= 2) & (df['ALT'].str.len() >= 2), ";Complex INDEL"),
        ((df['Review'] == "No Review") & (df['average_af'] >= 0.2), "High VAF"),
        ((df['Review'] != "No Review") & (df['average_af'] >= 0.2), ";High VAF"),
        ((df['Review'] == "No Review") & (df['Gene'].isin(gene_list)) & (df['VariantClass'].isin(missense_mutation)) & (df['homopolymerCase'] != ""), "Homopolymer Region"),
        ((df['Review'] != "No Review") & (df['Gene'].isin(gene_list)) & (df['VariantClass'].isin(missense_mutation)) & (df['homopolymerCase'] != ""), ";Homopolymer Region"),
        ((df['Review'] == "No Review") & (df['Gene'].isin(gene_list)) & (df['VariantClass'].isin(missense_mutation)) & (df['SIFT'].str.contains("deleterious")) & (df['PolyPhen'].str.contains("damaging")) & (df['putative_driver'] == 0) & (df['n.HGVSc'] >= 1), "B/B Missense Review"),
        ((df['Review'] != "No Review") & (df['Gene'].isin(gene_list)) & (df['VariantClass'].isin(missense_mutation)) & (df['SIFT'].str.contains("deleterious")) & (df['PolyPhen'].str.contains("damaging")) & (df['putative_driver'] == 0) & (df['n.HGVSc'] >= 1), ";B/B Missense Review"),
        ((df['Review'] == "No Review") & (df['Gene'].isin(gene_list)) & (df['VariantClass'].isin(missense_mutation)) & (df['SIFT'].str.contains("deleterious")) & (df['PolyPhen'].str.contains("damaging")) & (df['putative_driver'] == 0), "S/P Missense Review"),
        ((df['Review'] != "No Review") & (df['Gene'].isin(gene_list)) & (df['VariantClass'].isin(missense_mutation)) & (df['SIFT'].str.contains("deleterious")) & (df['PolyPhen'].str.contains("damaging")) & (df['putative_driver'] == 0), ";S/P Missense Review"),
        ((df['Review'] == "No Review") & (df['Gene'].isin(gene_list)) & (df['VariantClass'].isin(missense_mutation)) & ((df['nearBBLogic']) | (df['nearCosmicHemeLogic'])) & ((df['SIFT'].str.contains("deleterious")) | (df['PolyPhen'].str.contains("damaging"))) & (df['putative_driver'] == 0), "NHS Missense Review"),
        ((df['Review'] != "No Review") & (df['Gene'].isin(gene_list)) & (df['VariantClass'].isin(missense_mutation)) & ((df['nearBBLogic']) | (df['nearCosmicHemeLogic'])) & ((df['SIFT'].str.contains("deleterious")) | (df['PolyPhen'].str.contains("damaging"))) & (df['putative_driver'] == 0), ";NHS Missense Review"),
        ((df['Review'] == "No Review") & (df['Gene'].isin(gene_list)) & (df['VariantClass'].isin(missense_mutation)) & ((df['n.loci.vep'] - df['n.loci.truncating.vep']) >= 5) & ((df['SIFT'].str.contains("deleterious")) | (df['PolyPhen'].str.contains("damaging"))) & (df['putative_driver'] == 0), "HS Missense Review"),
        ((df['Review'] != "No Review") & (df['Gene'].isin(gene_list)) & (df['VariantClass'].isin(missense_mutation)) & ((df['n.loci.vep'] - df['n.loci.truncating.vep']) >= 5) & ((df['SIFT'].str.contains("deleterious")) | (df['PolyPhen'].str.contains("damaging"))) & (df['putative_driver'] == 0), ";HS Missense Review"),
        ((df['Review'] == "No Review") & (df['n_samples'] > 5) & (df['n.HGVSc'] < 25) & (df['CosmicCount'] < 50), "Recurrent"),
        ((df['Review'] != "No Review") & (df['n_samples'] > 5) & (df['n.HGVSc'] < 25) & (df['CosmicCount'] < 50), ";Recurrent"),
        ((df['Review'] == "No Review") & (df['Gene'].isin(unique_genes)), "Bick's Email"),
        ((df['Review'] != "No Review") & (df['Gene'].isin(unique_genes)), ";Bick's Email")
    ]
    for condition, label in conditions:
        df['Review'] = np.where(condition, label, df['Review'])

    # If something is ONLY in Review because it was recurrent, we can auto pass it if the recurrence wasn't that significant
    # 2.69745405 - chr4:105236829:C:T       <-- These numbers come from UKBB
    # 0.471336441 - chr22:19275752:G:A
    df['auto_pass_recurrence'] = np.where(
        (df['Review'].str.contains("Recurrent")) &
        (df['n.HGVSc'] == 0) &
        (df['heme_cosmic_count'] == 0) &
        ((df['ratio_to_BB'] <= 2.7) | (df['ratio_to_cosmic'] <= 0.47)),
        True,
        False
    )
    df['Review'] = np.where(df['auto_pass_recurrence'], df['Review'].str.replace("Recurrent", "Was Recurrent"), df['Review'])
    df = df.drop(columns=['auto_pass_recurrence', 'Protein_position_start', 'Protein_position_end'])

    if 'SpliceAI_pred_SYMBOL' in df.columns:
        df.loc[
            (df['Review'].str.contains("Splice Region Variant")) &
            (df['pd_reason'] == "Not PD") &
            ((df['SpliceAI_pred_DS_AG'] >= 0.8) |
            (df['SpliceAI_pred_DS_AL'] >= 0.8) |
            (df['SpliceAI_pred_DS_DG'] >= 0.8) |
            (df['SpliceAI_pred_DS_DL'] >= 0.8)),
            'pd_reason'
        ] = "Splice AI"
        df['putative_driver'] = np.where(df['pd_reason'] != "Not PD", 1, 0)

    review_df = df.loc[
        ((df['Review'].str.contains("Long INDEL")) & (df['putative_driver'] == 1)) |
        ((df['Review'].str.contains("Complex INDEL")) & (df['putative_driver'] == 1)) |
        ((df['Review'].str.contains("High VAF")) & (df['putative_driver'] == 1)) |
        ((df['Review'] == "Recurrent") & (df['putative_driver'] == 1)) |
        ((df['Review'].str.contains("Homopolymer Region")) & (df['putative_driver'] == 1)) |
        ((df['Review'].str.contains("Bick's Email")) & (df['putative_driver'] == 1)) |
        ((df['Review'].str.contains("Missense Review")) & (df['putative_driver'] == 0)) |
        ((df['Review'].str.contains("MW Review")) & (df['putative_driver'] == 0)) |
        ((df['Review'].str.contains("Splice Region Variant")) & (df['putative_driver'] == 0))
    ]

    pass_df = df.loc[
        ((df['putative_driver'] == 1) & (df['Review'] == "No Review")) |
        ((df['putative_driver'] == 1) & (df['Review'].str.contains("Missense Review"))) |
        ((df['putative_driver'] == 1) & (df['Review'] == "Was Recurrent")) |
        ((df['putative_driver'] == 1) & (df['Review'].str.contains("MW Review"))) |
        ((df['putative_driver'] == 1) & (df['Review'].str.contains("Splice Region Variant")))
    ]

    df['status'] = ""
    df['status'] = np.where(df['sample_key'].isin(review_df['sample_key']), "Review", df['status'])
    df['status'] = np.where(df['sample_key'].isin(pass_df['sample_key']), "Pass", df['status'])

    log.logit("Finished determining pathogenicity.")
    return review_df, pass_df, df


# CH Definition
# - ch_pd == 1
# - Pass gnomAD Filter
# - Save ASXL1 W646 Mutations
# - PASS Mutect + PASS Vardict
# - IF Mutect Filter == weak_evidence | strand_bias + HOTSPOT DEFINITION THEN PASS
# - Pass PoN2at2%
# - Pass Min VAF of 0.1%
# - max(mutect_vaf, vardict_vaf) >= 0.02
# - Mutect Min Alt >= 2 OR Vardict Min Alt >= 2
# - Mutect 1 AltCount on Fwd + Rev OR Vardict 1 AltCount on Fwd + Rev
# - Average VAF <= 0.25 OR ELSE B/B >= 5 OR ELSE CosmicCount >= 25
# - Save Nonsense Mutations in DTAP
# - Median VAF <= 0.35 if Average VAF > 0.25 (At least 2 Samples) <- For Tumors add the B/B and COSMIC
# - Calculate N Samples
# - PoN Edge Case of 0
def ch_to_df(mutect_connection, vardict_db, annotation_db, debug):
    log.logit(f"Grabbing CH Variants from Database...")
    mutect_connection.execute(f"ATTACH \'{vardict_db}\' as vardict_db")
    mutect_connection.execute(f"ATTACH \'{annotation_db}\' as annotation_db")
    sql = """
    CREATE OR REPLACE VIEW pd_filtered AS
    SELECT *
    FROM annotation_db.pd as p
    LEFT JOIN annotation_db.vep as vep
    ON p.variant_id = vep.variant_id
    WHERE (
        (max_gnomADe_AF_VEP < 0.005 OR max_gnomADe_AF_VEP is NULL) AND
        (max_gnomADg_AF_VEP < 0.005 OR max_gnomADg_AF_VEP is NULL) AND
        (max_pop_gnomAD_AF < 0.0005 OR max_pop_gnomAD_AF is NULL)
    ) OR (vep.key = 'chr20:32434638:A:AG' OR vep.key = 'chr20:32434638:A:AGG');

    CREATE OR REPLACE VIEW mutect_filtered AS
    SELECT mutect.*
    FROM mutect
    LEFT JOIN pd_filtered p
    ON mutect.variant_id = p.variant_id
    WHERE (
        (
            mutect_filter = '[PASS]' OR (
                (
                    mutect_filter = '[weak_evidence]' OR
                    mutect_filter = '[strand_bias]' OR
                    mutect_filter = '[weak_evidence, strand_bias]' OR
                    mutect_filter = '[strand_bias, weak_evidence]'
                ) AND
                (
                    \"n.loci.vep\" >= 5 OR
                    (CosmicCount >= 25 AND myeloid_cosmic_count >= 1) OR
                    CosmicCount >= 100 OR
                    heme_cosmic_count >= 10 OR
                    myeloid_cosmic_count >= 5 OR
                    (isTruncatingHotSpot == 1 AND (SYMBOL = 'DNMT3A' OR SYMBOL = 'TET2' OR SYMBOL = 'ASXL1' OR SYMBOL = 'PPM1D'))
                ) 
            )
        ) AND
            pon_2at2_percent is NULL AND
            format_af >= 0.001 AND
            mutect.variant_id IN (
                SELECT variant_id
                FROM pd_filtered
            )
        ) OR (mutect.key = 'chr20:32434638:A:AG' OR mutect.key = 'chr20:32434638:A:AGG');

    CREATE OR REPLACE VIEW vardict_filtered AS
    SELECT *
    FROM vardict_db.vardict
    WHERE (
        vardict_filter = '[PASS]' AND
        pon_2at2_percent is NULL AND
        format_af >= 0.001 AND
        variant_id IN (
            SELECT variant_id
            FROM pd_filtered
        )
    ) OR (key = 'chr20:32434638:A:AG' OR key = 'chr20:32434638:A:AGG');

    CREATE OR REPLACE VIEW n_samples_and_median_af AS
    SELECT m.variant_id, count(*) as n_samples, 
        MEDIAN((m.format_af + v.format_af)/2) AS median_af
    FROM mutect_filtered m
    INNER JOIN vardict_filtered v
    ON m.variant_id = v.variant_id AND m.sample_id = v.sample_id
    GROUP BY m.variant_id;

    CREATE OR REPLACE VIEW ch_pd AS
    SELECT m.sample_name, m.key, m.mutect_filter, m.info_mbq_ref, m.info_mbq_alt, m.info_mmq_ref, m.info_mmq_alt, m.format_af, m.format_dp, m.format_ref_count, m.format_alt_count, m.format_ref_f1r2, m.format_alt_f1r2, m.format_ref_f2r1, m.format_alt_f2r1, m.format_ref_fwd, m.format_ref_rev, m.format_alt_fwd, m.format_alt_rev, m.fisher_p_value, m.batch,
    v.vardict_filter, v.info_qual, v.info_mq, v.info_nm, v.format_ref_count, v.format_alt_count, v.format_dp, v.format_vd, v.format_af, v.format_ref_fwd, v.format_ref_rev, v.format_alt_fwd, v.format_alt_rev, v.fisher_p_value,
    (m.format_af + v.format_af)/2 as average_af, a.*
    FROM mutect_filtered m
    INNER JOIN vardict_filtered v
    ON m.variant_id = v.variant_id AND m.sample_id = v.sample_id
    LEFT JOIN pd_filtered a
    ON m.variant_id = a.variant_id
    WHERE 
        (GREATEST(m.format_af, v.format_af) >= 0.02) 
        AND
        (
            (m.format_alt_fwd >= 1 AND m.format_alt_rev >= 1) 
            OR 
            (v.format_alt_fwd >= 1 AND v.format_alt_rev >= 1)
        );

    SELECT c.*, s.n_samples, s.median_af
        FROM ch_pd c
        LEFT JOIN n_samples_and_median_af s
        ON c.variant_id = s.variant_id
        WHERE
            (
                average_af <= 0.25 OR \"n.HGVSc\" >= 5 OR CosmicCount >= 25 OR
                (
                    (SYMBOL = 'DNMT3A' OR SYMBOL = 'TET2' OR SYMBOL = 'ASXL1' OR SYMBOL = 'PPM1D') AND 
                    (VariantClass = 'frameshift_variant' OR VariantClass = 'stop_lost' OR VariantClass = 'stop_gained' OR VariantClass = 'transcript_ablation')
                )
            ) AND NOT (
                average_af >= 0.25 AND (median_af >= 0.35 AND n_samples > 1)
            )
    """
    #COPY () TO 'ch_pd.csv' (HEADER, DELIMITER ',');
    #AND fisher_p_value <= 1.260958e-09
    # NOT (
    #     (format_alt_fwd / (format_alt_fwd + format_alt_rev)) > 0.9 OR (format_alt_fwd / (format_alt_fwd + format_alt_rev)) < 0.1 AND
    #     (format_alt_rev / (format_alt_fwd + format_alt_rev)) > 0.9 OR (format_alt_rev / (format_alt_fwd + format_alt_rev))
    # )

    if debug: log.logit(f"Executing: {sql}")
    df = mutect_connection.execute(sql).df()
    mutect_connection.execute(f"DROP VIEW pd_filtered; DROP VIEW mutect_filtered; DROP VIEW vardict_filtered; DROP VIEW n_samples_and_median_af; DROP VIEW ch_pd")
    mutect_connection.execute(f"DETACH vardict_db")
    mutect_connection.execute(f"DETACH annotation_db")
    if debug: log.logit(f"SQL Complete")
    length = len(df)
    log.logit(f"{length} variants are identified to be CH mutations")
    return df

def dump_ch_variants(mutect_db, vardict_db, annotation_db, debug):
    connection = db.duckdb_connect_rw(mutect_db, False)
    df = ch_to_df(connection, vardict_db, annotation_db, debug)
    total_sample = connection.execute(f"SELECT COUNT(DISTINCT sample_id) FROM mutect").fetchone()[0]
    review_df, pass_df, df = determine_pathogenicity(df, total_sample, debug)
    length = len(df)
    print(df)
    df.to_csv(f"ch_pd.all.csv", index=False, mode='w')
    review_df.to_csv(f"ch_pd.review.csv", index=False, mode='w')
    pass_df.to_csv(f"ch_pd.pass.csv", index=False, mode='w')
    log.logit(f"{length} variants are putative drivers")
    return df
