import decimal
import scipy.stats as stats
from scipy.stats import fisher_exact
import numpy as np

def pvalue(PoN_RefDepth, PoN_AltDepth, RefDepth, AltDepth):
    if PoN_RefDepth is None or PoN_AltDepth is None:
        return None
    else:
        tbl = [[PoN_RefDepth, PoN_AltDepth], [RefDepth, AltDepth]]
        return fisher_exact(tbl, alternative='two-sided').pvalue

def pvalue_variant_sample(variant_id, sample_id, ref_fwd, ref_rev, alt_fwd, alt_rev, PoN_RefDepth, PoN_AltDepth):
    tbl = [[PoN_RefDepth, PoN_AltDepth], [ref_fwd + ref_rev, alt_fwd + alt_rev]]
    res = fisher_exact(tbl, alternative='two-sided')
    return res.pvalue, variant_id, sample_id

# Fisher Test using fisher package (much faster)
def pvalue_df(df):
    df['rd'] = df['format_ref_fwd'] + df['format_ref_rev']
    df['ad'] = df['format_alt_fwd'] + df['format_alt_rev']
    c = df[['PoN_RefDepth','PoN_AltDepth','rd','ad']].to_numpy(dtype='uint')
    from fisher import pvalue_npy
    _, _, twosided = pvalue_npy(c[:, 0], c[:, 1], c[:, 2], c[:, 3])
    df['pvalue'] = twosided
    return(df)
