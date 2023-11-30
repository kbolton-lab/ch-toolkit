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
    # Checking for conditions
    # 1. When PoN_AltDepth is 0
    # 2. When PoN VAF >= Variant VAF
    df.loc[(df['PoN_AltDepth'] == 0) & (df['PoN_RefDepth'] > 0), 'pvalue'] = 0
    df.loc[(df['PoN_AltDepth'] / (df['PoN_AltDepth'] + df['PoN_RefDepth']) >= df['ad'] / (df['ad'] + df['rd'])), 'pvalue'] = 1
    df.loc[((df['PoN_AltDepth'] == 0) & (df['PoN_AltDepth'] != 0)) & ((df['rd'] == 0) & df['ad'] != 0), 'pvalue'] = 1
    return(df)