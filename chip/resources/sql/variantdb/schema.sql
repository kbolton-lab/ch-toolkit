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

CREATE TABLE IF NOT EXISTS fp_filter_types
(
    fp_filter_type_id     integer NOT NULL PRIMARY KEY,
    tag                   varchar(10) NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS variants
(
    variant_id             integer PRIMARY KEY,
    chrom_id               integer NOT NULL REFEFERENCES chromosomes(chrom_id),
    pos                    integer NOT NULL,
    ref                    varchar(255) NOT NULL,
    alt                    varchar(255) NOT NULL,
    qc_pass                boolean CHECK (qc_pass in (0,1)),
    UNIQUE( chrom_id, pos, ref, alt )
);

CREATE TABLE IF NOT EXISTS fp_filter
(
    sample_id              integer NOT NULL REFERENCES samples(sample_id),
    variant_id             integer NOT NULL REFERENCES variants(variant_id),
    fp_filter_type_id      integer REFERENCES fp_filter_types(fp_filter_type_id),
    PRIMARY KEY( sample_id, variant_id, fp_filter_type_id )
);

CREATE TABLE IF NOT EXISTS mutect
(
    sample_id              integer NOT NULL REFERENCES samples(sample_id),
    variant_id             integer NOT NULL REFERENCES variants(variant_id),
    version                varchar(10),
    filter_foo             varchar(10),
    filter_bar             varchar(16),
    info_x                 varchar(16),
    info_y                 varchar(16),
    format_ref_count       integer,
    format_alt_count       integer,
    fisher_p_value         decimal(10,5), /* indicates that the field should be used to store a value up to ten digits in length, with up to five digits before the decimal point and up to five digits after the decimal point. */
    PRIMARY KEY( sample_id, variant_id )
);

CREATE TABLE IF NOT EXISTS mutect
(
    sample_id              integer NOT NULL REFERENCES samples(sample_id),
    variant_id             integer NOT NULL REFERENCES variants(variant_id),
    version                varchar(10),
    filter_foo             varchar(10),
    filter_bar             varchar(16),
    info_x                 varchar(16),
    info_y                 varchar(16),
    format_ref_count       integer,
    format_alt_count       integer,
    fisher_p_value         decimal(10,5),
    PRIMARY KEY( sample_id, variant_id )
);

CREATE TABLE IF NOT EXISTS vardict
(
    sample_id              integer NOT NULL REFERENCES samples(sample_id),
    variant_id             integer NOT NULL REFERENCES variants(variant_id),
    version                varchar(10),
    filter_foo             varchar(10),
    filter_bar             varchar(16),
    info_x                 varchar(16),
    info_y                 varchar(16),
    format_ref_count       integer,
    format_alt_count       integer,
    fisher_p_value         decimal(10,5),
    PRIMARY KEY( sample_id, variant_id )
);

CREATE TABLE IF NOT EXISTS lofreq
(
    sample_id              integer NOT NULL REFERENCES samples(sample_id),
    variant_id             integer NOT NULL REFERENCES variants(variant_id),
    version                varchar(10),
    filter_foo             varchar(10),
    filter_bar             varchar(16),
    info_x                 varchar(16),
    info_y                 varchar(16),
    format_ref_count       integer,
    format_alt_count       integer,
    fisher_p_value         decimal(10,5),
    PRIMARY KEY( sample_id, variant_id )
);

CREATE TABLE IF NOT EXISTS pindel
(
    sample_id              integer NOT NULL REFERENCES samples(sample_id),
    variant_id             integer NOT NULL REFERENCES variants(variant_id),
    version                varchar(10),
    filter_foo             varchar(10),
    filter_bar             varchar(16),
    info_x                 varchar(16),
    info_y                 varchar(16),
    format_ref_count       integer,
    format_alt_count       integer,
    fisher_p_value         decimal(10,5),
    PRIMARY KEY( sample_id, variant_id )
);
