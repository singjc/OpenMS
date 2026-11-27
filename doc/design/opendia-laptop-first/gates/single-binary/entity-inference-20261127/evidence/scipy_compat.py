"""Value-preserving SciPy 1.18 rankdata dtype compatibility for pinned PyProphet.

SciPy 1.18.1 returns float64 for method='min'. Pinned PyProphet uses
that array directly as a NumPy index in stats.lfdr. Restore only integer rank
storage for method='min', never change values, tie rules or estimator options.
The upstream repository remains read-only and byte-verified against its pin.
"""
import numpy as np
import scipy
import scipy.stats

_original_rankdata = scipy.stats.rankdata


def rankdata_compat(a, method="average", *args, **kwargs):
    ranks = _original_rankdata(a, method, *args, **kwargs)
    if method == "min" and np.all(np.isfinite(ranks)):
        integers = ranks.astype(np.int64)
        np.testing.assert_array_equal(ranks, integers)
        return integers
    return ranks


def enable():
    scipy.stats.rankdata = rankdata_compat
    return {"scipy": scipy.__version__, "shim": "rankdata(method=min) integer dtype only, identical rank values"}
