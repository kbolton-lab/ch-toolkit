CREATE TABLE IF NOT EXISTS samples
(
    sample_id          integer NOT NULL PRIMARY KEY,
    sample_name        text NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS chromosomes
(
    chrom_id     integer PRIMARY KEY,
    chrom        varchar(5) NOT NULL UNIQUE
);


CREATE TABLE IF NOT EXISTS variants
(
    variant_id             integer PRIMARY KEY,
    chrom_id               integer NOT NULL REFEFERENCES chromosomes(chrom_id),
    pos                    integer NOT NULL,
    ref                    varchar(255) NOT NULL,
    alt                    varchar(255) NOT NULL,
    pon_refdepth           integer,
    pon_altdepth           integer,
    pon_2at2_percent       boolean CHECK (pon_2at2_percent in (0,1)),
    pon_nat2_percent       integer,
    pon_max_vaf            decimal(10,5),
    vep                    varchar(10000),
    qc_pass                boolean CHECK (qc_pass in (0,1)),
    UNIQUE( chrom_id, pos, ref, alt )
);

/* VEP EXAMPLE
CSQ=-|frameshift_variant|HIGH|DNMT3A|ENSG00000119772|Transcript|ENST00000321117.10|protein_coding|17/23||ENST00000321117.10:c.1954del|ENSP00000324375.5:p.Asp652ThrfsTer53|2231|1954|652|D/X|Gac/ac|||-1||deletion|HGNC|HGNC:2978|YES|NM_022552.5|1|P3|CCDS33157.1|ENSP00000324375|Q9Y6K1.187||UPI000000DA70||Ensembl||C|C||1|||PDB-ENSP_mappings:2qrv.A&PDB-ENSP_mappings:2qrv.D&PDB-ENSP_mappings:2qrv.E&PDB-ENSP_mappings:2qrv.H&Gene3D:3.40.50.150&PDB-ENSP_mappings:4u7p.A&PDB-ENSP_mappings:4u7t.A&PDB-ENSP_mappings:4u7t.C&PDB-ENSP_mappings:5yx2.A&PDB-ENSP_mappings:5yx2.D&PDB-ENSP_mappings:6brr.A&PDB-ENSP_mappings:6brr.D&PDB-ENSP_mappings:6f57.A&PDB-ENSP_mappings:6f57.D&PDB-ENSP_mappings:6pa7.K&PDB-ENSP_mappings:6pa7.P&PDB-ENSP_mappings:6w89.A&PDB-ENSP_mappings:6w89.D&PDB-ENSP_mappings:6w89.G&PDB-ENSP_mappings:6w89.J&PDB-ENSP_mappings:6w8b.A&PDB-ENSP_mappings:6w8b.D&PDB-ENSP_mappings:6w8b.H&PDB-ENSP_mappings:6w8b.K&PDB-ENSP_mappings:6w8d.A&PDB-ENSP_mappings:6w8d.D&PDB-ENSP_mappings:6w8j.A&PDB-ENSP_mappings:6w8j.D&Pfam:PF00145&PROSITE_profiles:PS51679&PANTHER:PTHR23068&PANTHER:PTHR23068:SF10&Superfamily:SSF53335||||||||||||||||||||||||||||||||MPAMPSSGPGDTSSSAAEREEDRKDGEEQEEPRGKEERQEPSTTARKVGRPGRKRKHPPVESGDTPKDPAVISKSPSMAQDSGASELLPNGDLEKRSEPQPEEGSPAGGQKGGAPAEGEGAAETLPEASRAVENGCCTPKEGRGAPAEAGKEQKETNIESMKMEGSRGRLRGGLGWESSLRQRPMPRLTFQAGDPYYISKRKRDEWLARWKREAEKKAKVIAGMNAVEENQGPGESQKVEEASPPAVQQPTDPASPTVATTPEPVGSDAGDKNATKAGDDEPEYEDGRGFGIGELVWGKLRGFSWWPGRIVSWWMTGRSRAAEGTRWVMWFGDGKFSVVCVEKLMPLSSFCSAFHQATYNKQPMYRKAIYEVLQVASSRAGKLFPVCHDSDESDTAKAVEVQNKPMIEWALGGFQPSGPKGLEPPEEEKNPYKEVYTDMWVEPEAAAYAPPPPAKKPRKSTAEKPKVKEIIDERTRERLVYEVRQKCRNIEDICISCGSLNVTLEHPLFVGGMCQNCKNCFLECAYQYDDDGYQSYCTICCGGREVLMCGNNNCCRCFCVECVDLLVGPGAAQAAIKEDPWNCYMCGHKGTYGLLRRREDWPSRLQMFFANNHDQEFDPPKVYPPVPAEKRKPIRVLSLFDGIATGLLVLKTWAFRWTATLPRRCVRTPSRWAWCGTRGRSCTSGTSAASHRSISRSGAHSIW|MPAMPSSGPGDTSSSAAEREEDRKDGEEQEEPRGKEERQEPSTTARKVGRPGRKRKHPPVESGDTPKDPAVISKSPSMAQDSGASELLPNGDLEKRSEPQPEEGSPAGGQKGGAPAEGEGAAETLPEASRAVENGCCTPKEGRGAPAEAGKEQKETNIESMKMEGSRGRLRGGLGWESSLRQRPMPRLTFQAGDPYYISKRKRDEWLARWKREAEKKAKVIAGMNAVEENQGPGESQKVEEASPPAVQQPTDPASPTVATTPEPVGSDAGDKNATKAGDDEPEYEDGRGFGIGELVWGKLRGFSWWPGRIVSWWMTGRSRAAEGTRWVMWFGDGKFSVVCVEKLMPLSSFCSAFHQATYNKQPMYRKAIYEVLQVASSRAGKLFPVCHDSDESDTAKAVEVQNKPMIEWALGGFQPSGPKGLEPPEEEKNPYKEVYTDMWVEPEAAAYAPPPPAKKPRKSTAEKPKVKEIIDERTRERLVYEVRQKCRNIEDICISCGSLNVTLEHPLFVGGMCQNCKNCFLECAYQYDDDGYQSYCTICCGGREVLMCGNNNCCRCFCVECVDLLVGPGAAQAAIKEDPWNCYMCGHKGTYGLLRRREDWPSRLQMFFANNHDQEFDPPKVYPPVPAEKRKPIRVLSLFDGIATGLLVLKDLGIQVDRYIASEVCEDSITVGMVRHQGKIMYVGDVRSVTQKHIQEWGPFDLVIGGSPCNDLSIVNPARKGLYEGTGRLFFEFYRLLHDARPKEGDDRPFFWLFENVVAMGVSDKRDISRFLESNPVMIDAKEVSAAHRARYFWGNLPGMNRPLASTVNDKLELQECLEHGRIAKFSKVRTITTRSNSIKQGKDQHFPVFMNEKEDILWCTEMERVFGFPVHYTDVSNMSRLARQRLLGRSWSVPVIRHLFAPLKEYFACV|-1|-15|-38|-14|0.00|0.12|0.01|0.00|DNMT3A||||||||||||||||||||||||||||||
*/

CREATE TABLE IF NOT EXISTS fp_filter_types
(
    fp_filter_type_id     integer NOT NULL PRIMARY KEY,
    tag                   varchar(10) NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS fp_filter
(
    sample_id              integer NOT NULL REFERENCES samples(sample_id),
    variant_id             integer NOT NULL REFERENCES variants(variant_id),
    fp_filter_type_id      integer REFERENCES fp_filter_types(fp_filter_type_id),
    PRIMARY KEY( sample_id, variant_id, fp_filter_type_id )
);

CREATE TABLE IF NOT EXISTS annotate_pd
(
    sample_id                       integer NOT NULL REFERENCES samples(sample_id),
    variant_id                      integer NOT NULL REFERENCES variants(variant_id),
    max_gnomAD_AF_VEP               decimal(2,10),
    max_gnomADe_AF_VEP              decimal(2,10),
    max_gnomADg_AF_VEP              decimal(2,10),
    max_gnomAD_AF                   boolean CHECK (max_gnomAD_AF in (0,1)),
    context_5                       varchar(10),
    context_3                       varchar(10),
    context_10                      varchar(10),
    dust_score_5                    decimal(2,10),
    dust_score_3                    decimal(2,10),
    dust_score_10                   decimal(2,10),
    dust_score                      decimal(2,10),
    case_NXXX                       boolean CHECK (case_NXXX in (0,1)),
    case_XNXX                       boolean CHECK (case_XNXX in (0,1)),
    case_XXNX                       boolean CHECK (case_XXNX in (0,1)),
    case_XXXN                       boolean CHECK (case_XXXN in (0,1)),
    case_NNXX                       boolean CHECK (case_NNXX in (0,1)),
    case_XNNX                       boolean CHECK (case_XNNX in (0,1)),
    case_XXNN                       boolean CHECK (case_XXNN in (0,1)),
    deletion                        varchar(10000),
    deletion_plus_context_3         varchar(10000),
    case_deletion_in_homopolymer    boolean CHECK (case_deletion_in_homopolymer in (0,1))
    pass_homopolymer_filter         boolean CHECK (pass_homopolymer_filter in (0,1))
    pass_strand_bias_mutect         boolean CHECK (pass_strand_bias_mutect in (0,1)),
    pass_strand_bias_vardict        boolean CHECK (pass_strand_bias_vardict in (0,1)),
    pass_strand_bias_lofreq         boolean CHECK (pass_strand_bias_lofreq in (0,1)),
    AAchange                        varchar(100),
    gene_loci_p                     varchar(100),
    gene_loci_c                     varchar(100),
    gene_loci_vep                   varchar(100),       /* ASXL1_Y591_R678delinsFQGCSGLPWAQ */
    gene_aachange                   varchar(100),
    gene_cDNAchange                 varchar(100),
    n_loci_vep                      integer,
    source_totals_loci              varchar(30),
    n_loci_truncating_vep           integer,
    source_totals_loci_truncating   varchar(30),
    n_HGVSp                         integer,
    source_totals_p                 varchar(30),
    n_HGVSc                         integer,
    source_totals_c                 varchar(30),
    COSMIC_ID                       varchar(12),
    CosmicCount                     integer,
    heme_cosmic_count               integer,
    myeloid_cosmic_count            integer,
    oncoKB                          varchar(30),
    isOncogenic                     boolean CHECK (isOncogenic in (0,1))
    isTSG                           boolean CHECK (isTSG in (0,1))
    isTruncatingHotSpot             boolean CHECK (isTruncatingHotSpot in (0,1))
    ch_my_pd                        integer,
    ch_pd                           integer,
    ch_pd2                          integer,
    VariantClass                    varchar(20),
    Gene                            varchar(10),
    WHY_CH                          varchar(100),           /* "{""ch_my_pd"":[""PD_table;""],""ch_pd2"":[""CH_pd=1;""],""ch_pd"":[""OncoKB; ch_my_pd>=1""]}" */
    PRIMARY KEY( sample_id, variant_id )
);

CREATE TABLE IF NOT EXISTS mutect_filter_types
(
    mutect_filter_type_id integer NOT NULL PRIMARY KEY,
    tag                   varchar(10) NOT NULL UNIQUE
);


CREATE TABLE IF NOT EXISTS mutect
(
    sample_id                           integer NOT NULL REFERENCES samples(sample_id),
    variant_id                          integer NOT NULL REFERENCES variants(variant_id),
    version                             varchar(10),
    mutect_filter_type_id               integer REFERENCES mutect_filter_types(mutect_filter_type_id),
    info_as_filterstatus                varchar(20),     /* e.g. AS_FilterStatus=strand_bias|strand_bias */
    info_as_sb_table                    varchar(20),     /* e.g. AS_SB_TABLE=5379,6627|1,27|1,174 */
    info_as_uniq_alt_read_count         integer,
    info_contq                          decimal(10,1),
    info_dp                             integer,
    info_ecnt                           integer,
    info_mbq_ref                        integer,
    info_mbq_alt                        integer,
    info_mfrl_ref                       integer,
    info_mfrl_alt                       integer,
    info_mmq_ref                        integer,
    info_mmq_alt                        integer,
    info_mpos                           integer,
    info_nalod                          decimal(10,1),
    info_ncount                         integer,
    info_nlod                           decimal(10,1),
    info_ocm                            integer,
    info_pon                            varchar(10),     /* 0 - This denotes that this value needs to be in double quotations*/
    info_popaf                          decimal(10,1),
    info_roq                            decimal(10,1),
    info_rpa_ref                        integer,
    info_rpa_alt                        integer,
    info_ru                             varchar(10),
    info_seqq                           integer,
    info_str                            varchar(10),    /* 0 */
    info_strandq                        integer,
    info_strq                           integer,
    info_tlod                           decimal(10,1),
    format_ref_count                    integer,
    format_alt_count                    integer,
    format_af                           decimal(10,5),
    format_dp                           integer,
    format_ref_f1r2                     integer,
    format_alt_f1r2                     integer,
    format_ref_f2r1                     integer,
    format_alt_f2r1                     integer,
    format_gq                           integer,
    format_gt                           varchar(3),     /* e.g. 0/1 */
    format_pgt                          varchar(5),     /* Not sure if relevant */
    format_pid                          varchar(5),     /* Not sure if relevant */
    format_pl                           integer,        /* G - Genotype */
    format_ps                           integer,
    format_sb                           varchar(20),    /* e.g. 14922,9228,61,28 */
    fisher_p_value                      decimal(10,5), /* indicates that the field should be used to store a value up to ten digits in length, with up to five digits before the decimal point and up to five digits after the decimal point. */
    info_old_multiallelic               varchar(100),  /* If this variant was multiallelic, saves original encoding */
    PRIMARY KEY( sample_id, variant_id )
);

CREATE TABLE IF NOT EXISTS vardict_filter_types
(
    vardict_filter_type_id  integer NOT NULL PRIMARY KEY,
    tag                     varchar(10) NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS vardict
(
    sample_id               integer NOT NULL REFERENCES samples(sample_id),
    variant_id              integer NOT NULL REFERENCES variants(variant_id),
    version                 varchar(10),
    vardict_filter_type_id  integer REFERENCES vardict_filter_types(vardict_filter_type_id),
    info_type               varchar(5),
    info_dp                 integer,
    info_end                integer,
    info_vd                 integer,
    info_af                 integer,
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
    info_gdamp              integer,
    info_tlamp              integer,
    info_ncamp              integer,
    info_ampflag            integer,
    info_hicnt              integer,
    info_hicov              integer,
    info_splitread          integer,
    info_spanpair           integer,
    info_svtype             varchar(3),
    info_svlen              integer,
    info_duprate            decimal(10,5),
    format_ref_count        integer,
    format_alt_count        integer,
    format_gt               varchar(3),
    format_dp               integer,
    format_cd               integer,
    format_af               decimal(10,5),
    format_ref_rd           integer,
    format_alt_rd           integer,
    format_ref_ald          integer,
    format_alt_ald          integer,
    fisher_p_value          decimal(10,5),
    info_old_multiallelic   varchar(100),  /* If this variant was multiallelic, saves original encoding */
    PRIMARY KEY( sample_id, variant_id )
);

CREATE TABLE IF NOT EXISTS lofreq_filter_types
(
    lofreq_filter_type_id  integer NOT NULL PRIMARY KEY,
    tag                     varchar(10) NOT NULL UNIQUE
);


CREATE TABLE IF NOT EXISTS lofreq
(
    sample_id              integer NOT NULL REFERENCES samples(sample_id),
    variant_id             integer NOT NULL REFERENCES variants(variant_id),
    version                varchar(10),
    lofreq_filter_type_id  integer REFERENCES lofreq_filter_types(lofreq_filter_type_id),
    info_dp                integer,
    info_af                decimal(10,5),
    info_sb                integer,
    info_dp_reffwd         integer,
    info_dp_refrev         integer,
    info_dp_altfwd         integer,
    info_dp_altrev         integer,
    info_indel             boolean,         /* If it is an INDEL, this flag just exists in the INFO column */
    info_consvar           boolean,         /* If it is an INDEL, this flag just exists in the INFO column */
    info_hrun              integer,
    format_ref_count       integer,
    format_alt_count       integer,
    format_gt              varchar(3),
    format_dp              integer,
    format_af              decimal(10,5),
    format_sb              integer,
    fisher_p_value         decimal(10,5),
    info_old_multiallelic  varchar(100),  /* If this variant was multiallelic, saves original encoding */
    PRIMARY KEY( sample_id, variant_id )
);

CREATE TABLE IF NOT EXISTS pindel
(
    sample_id              integer NOT NULL REFERENCES samples(sample_id),
    variant_id             integer NOT NULL REFERENCES variants(variant_id),
    version                varchar(10),
    info_homlen            integer,
    info_pf                integer,
    info_homseq            varchar(10),
    info_svlen             integer,
    info_svtype            varchar(3),
    info_ntlen             integer,
    format_gt              varchar(3),
    format_ref_count       integer,
    format_alt_count       integer,
    fisher_p_value         decimal(10,5),
    PRIMARY KEY( sample_id, variant_id )
);
