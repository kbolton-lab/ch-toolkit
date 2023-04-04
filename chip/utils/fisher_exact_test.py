import decimal
from scipy.stats import fisher_exact

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