import decimal
import numpy as np
from scipy.stats import fisher_exact
from numba import njit

# Fisher Test using fisher package (much faster)
def pvalue_df(df):
    df['rd'] = df['format_ref_fwd'] + df['format_ref_rev']
    df['ad'] = df['format_alt_fwd'] + df['format_alt_rev']

    # This method is faster, but has issues with accuracy with low p-values
    #c = df[['PoN_RefDepth','PoN_AltDepth','rd','ad']].to_numpy(dtype='uint')
    #from fisher import pvalue_npy
    #_, _, twosided = pvalue_npy(c[:, 0], c[:, 1], c[:, 2], c[:, 3])
    #df['pvalue'] = twosided
    #df['pvalue'] = [fisher_exact([[r[0], r[1]], [r[2], r[3]]])[1] for r in  df[['PoN_RefDepth','PoN_AltDepth','rd','ad']].values]

    maxSum = df[['PoN_RefDepth','PoN_AltDepth','rd','ad']].sum(axis = 1).max()
    logFac = unitLogFac(maxSum)
    df['pvalue'] = [ufet(r[0], r[1], r[2], r[3], logFac) for r in  df[['PoN_RefDepth','PoN_AltDepth','rd','ad']].values]

    # Checking for conditions
    # 1. When PoN_AltDepth is 0
    # 2. When PoN VAF >= Variant VAF
    df.loc[(df['PoN_AltDepth'] == 0) & (df['PoN_RefDepth'] > 0), 'pvalue'] = 0
    df.loc[(df['PoN_AltDepth'] / (df['PoN_AltDepth'] + df['PoN_RefDepth']) >= df['ad'] / (df['ad'] + df['rd'])), 'pvalue'] = 1
    df.loc[((df['PoN_AltDepth'] == 0) & (df['PoN_AltDepth'] != 0)) & ((df['rd'] == 0) & df['ad'] != 0), 'pvalue'] = 1
    return(df)

@njit
def ulogHypergeometricProb(logFacs, a, b, c, d):
    return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d] \
        - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d]

@njit
def unitLogFac(m):
    res = np.arange(0, m+1)  # np.array(range(1, m+1))
    res = np.log(res)
    res[0] = 0
    for i in np.arange(1, len(res)):
        res[i] += res[i-1]
    return res

@njit
def ufet(a, b, c, d, logFacs):
    n = a + b + c + d
    f = 10_000_000
    logpCutoff = round(ulogHypergeometricProb(logFacs, a, b, c, d) * f) / f
    pFraction = 0
    for x in range(n+1):
        if a+b-x >= 0 and a+c-x >= 0 and d-a+x >= 0:
            l = round(ulogHypergeometricProb(
                logFacs, x, a+b-x, a+c-x, d-a+x) * f) / f
            if l <= logpCutoff:
                pFraction += np.exp(l - logpCutoff)

    logpValue = logpCutoff + np.log(pFraction)
    pval = np.exp(logpValue)

    return pval