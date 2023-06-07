import decimal
import numpy as np

# Fisher Test using fisher package (much faster)
def pvalue_df(df):
    df['rd'] = df['format_ref_fwd'] + df['format_ref_rev']
    df['ad'] = df['format_alt_fwd'] + df['format_alt_rev']
    c = df[['PoN_RefDepth','PoN_AltDepth','rd','ad']].to_numpy(dtype='uint')
    from fisher import pvalue_npy
    _, _, twosided = pvalue_npy(c[:, 0], c[:, 1], c[:, 2], c[:, 3])
    df['pvalue'] = twosided
    return(df)
