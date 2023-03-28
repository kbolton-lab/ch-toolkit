"""Calculate the factorial of a number."""
def factorial(n):
    if n == 0: return 1
    else: return n * factorial(n-1)

"""Calculate the binomial coefficient n choose k."""
def choose(n, k):
    return factorial(n) / (factorial(k) * factorial(n-k))

def pvalue(PoN_RefDepth, PoN_AltDepth, RefDepth, AltDepth):
    if PoN_RefDepth is None or PoN_AltDepth is None:
        return 'NULL'
    else:
        tbl = [[PoN_RefDepth, PoN_AltDepth], [RefDepth, AltDepth]]
        return choose(sum(tbl[0]), tbl[0][0]) * choose(sum(tbl[1]), tbl[1][0]) / choose(sum(map(sum, tbl)), sum(tbl[0]))
