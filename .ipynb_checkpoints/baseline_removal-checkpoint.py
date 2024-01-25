import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve

# lam value should be between 10e2 and 10e9
# p value should be between 0.1 and 0.001
def baseline_als(y, lam=10e2, p=0.1, niter=10):
    """
    :params y intensity
    :params lam value should be between 10e2 and 10e9 - 10e2 by default
    :params p value should be between 0.1 and 0.001 - 0.1 by default
    :params niter = 10
    """
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z