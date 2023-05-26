import sys
import duckdb
import pandas as pd
import chip.utils.logger as log
import chip.utils.database as db

# CH Definition
# - Pass Mutect & Pass Vardict                          DONE
# - Pass PON, PON2,                                     DONE - Assume Fisher
# - max(mutect_vaf, vardict_vaf) >= 0.02
# - and gnomAD                                          DONE
# - Mutect Strand Alt >= 1 OR Vardict Strand Alt >= 1   DONE
# - Mutect Min Alt >= 2 OR Vardict Min Alt >= 2         DONE
# - Calculate MedianVAF

def dump_ch_variants(mutect_db, vardict_db, annotation_db, debug):
    connection = db.duckdb_connect_rw(mutect_db, False)
    connection.execute(f"ATTACH \'{vardict_db}\' as vardict_db")
    connection.execute(f"ATTACH \'{annotation_db}\' as annotation_db")
    sql = """
    CREATE OR REPLACE VIEW pd_filtered AS
    SELECT *
    FROM annotation_db.pd as p
    LEFT JOIN annotation_db.vep as vep
    ON p.variant_id = vep.variant_id
    WHERE (
        (max_gnomAD_AF_VEP < 0.005 OR max_gnomAD_AF_VEP is NULL) AND
        (max_gnomADe_AF_VEP < 0.005 OR max_gnomADe_AF_VEP is NULL) AND
        (max_gnomADg_AF_VEP < 0.005 OR max_gnomADg_AF_VEP is NULL) AND
        (max_pop_gnomAD_AF < 0.0005 OR max_pop_gnomAD_AF is NULL)
    )
    AND
    (
        \"n.loci.vep\" >= 5 OR
        (CosmicCount >= 25 AND myeloid_cosmic_count >= 1) OR
        CosmicCount >= 100 OR
        heme_cosmic_count >= 10 OR
        myeloid_cosmic_count >= 5 OR
        isTruncatingHotSpot == 1
    );

    CREATE OR REPLACE VIEW mutect_filtered AS
    SELECT *
    FROM mutect
    WHERE (
        (
            mutect_filter = '[PASS]' OR
            mutect_filter = '[weak_evidence]' OR
            mutect_filter = '[strand_bias]' OR
            mutect_filter = '[weak_evidence, strand_bias]' OR
            mutect_filter = '[strand_bias, weak_evidence]'
        ) AND
            pon_2at2_percent is NULL AND
            format_af >= 0.001 AND
            variant_id IN (
                SELECT variant_id
                FROM pd_filtered
            )
        ) OR (key = 'chr20:32434638:A:AG' OR key = 'chr20:32434638:A:AGG');

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

    SELECT *
        FROM mutect_filtered m
        INNER JOIN vardict_filtered v
        ON m.variant_id = v.variant_id AND m.sample_id = v.sample_id
        LEFT JOIN pd_filtered a
        ON m.variant_id = a.variant_id
        WHERE (m.format_alt_fwd >= 1 AND m.format_alt_rev >= 1) OR (v.format_alt_fwd >= 1 AND v.format_alt_rev >= 1);
    """

    #AND fisher_p_value <= 1.260958e-09
    #
    if debug: log.logit(f"Executing: {sql}")
    df = connection.execute(sql).df()
    connection.execute(f"DROP VIEW pd_filtered; DROP VIEW mutect_filtered; DROP VIEW vardict_filtered;")
    connection.execute(f"DETACH vardict_db")
    connection.execute(f"DETACH annotation_db")
    if debug: log.logit(f"SQL Complete")
    length = len(df)
    print(df)
    log.logit(f"{length} variants are putative drivers")
    return df

#SELECT max_gnomAD_AF_VEP, max_gnomADe_AF_VEP, max_gnomADg_AF_VEP, max_pop_gnomAD_AF, "n.loci.vep" FROM a.vep WHERE key = 'chr22:28710060:C:T';
#SELECT CosmicCount, myeloid_cosmic_count, heme_cosmic_count, isTruncatingHotSpot FROM a.pd WHERE key = 'chr22:28710060:C:T';
#SELECT mutect_filter, pon_2at2_percent, format_af, format_alt_fwd, format_alt_rev FROM mutect WHERE key = 'chr22:28710060:C:T';
#SELECT vardict_filter, pon_2at2_percent, format_af, format_alt_fwd, format_alt_rev FROM v.vardict WHERE key = 'chr22:28710060:C:T';
