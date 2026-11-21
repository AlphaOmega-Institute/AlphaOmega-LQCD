"""Charge of a chain of maps and its cohomology (from the BRST working notes, AlphaOmega-Theory,
manuscripts/notes/brst/test_brst_families.py). Q = sum_k P_{G_k} J with J = (-1)^degree; its reduced
space is the cohomology of the chain, with Krein sign alternating by degree."""
import numpy as np
from scipy.linalg import null_space, orth

TOL = 1e-9


def polar(d):
    U, s, Vh = np.linalg.svd(d)
    r = int((s > TOL).sum())
    return U[:, :r] @ Vh[:r], Vh[:r].conj().T @ Vh[:r], U[:, :r] @ U[:, :r].conj().T


def chain_charge(dims, maps, degrees):
    """Q = sum_k P_{G_k} J for the chain C_0 -> C_1 -> ... with J = (-1)^degree, and the ghost
    number N = (1/2) sum_k R_k built from the polar reflections of the maps."""
    off = np.cumsum([0] + list(dims))
    n = off[-1]
    J = np.diag(np.concatenate([[(-1.0) ** degrees[i]] * dims[i] for i in range(len(dims))]))
    PG = np.zeros((n, n), complex)
    R = np.zeros((n, n), complex)
    for i, d in enumerate(maps):
        V, PS, PD = polar(d)
        a, b = slice(off[i], off[i + 1]), slice(off[i + 1], off[i + 2])
        PG[a, a] += PS / 2
        PG[a, b] += V.conj().T / 2
        PG[b, a] += V / 2
        PG[b, b] += PD / 2
        R[a, b] += V.conj().T
        R[b, a] += V
    return PG @ J, J, R / 2


def cohomology(Q, J):
    K = null_space(Q, rcond=TOL)
    Rg = orth(Q, rcond=TOL)
    H = K @ null_space(Rg.conj().T @ K, rcond=TOL) if Rg.size else K
    w = np.linalg.eigvalsh(H.conj().T @ J @ H) if H.size else np.array([])
    return H, int((w > TOL).sum()), int((w < -TOL).sum())


def chain_cohomology_dims(dims, maps):
    out = []
    for i in range(len(dims)):
        ker = dims[i] - (np.linalg.matrix_rank(maps[i], tol=TOL) if i < len(maps) else 0)
        im = np.linalg.matrix_rank(maps[i - 1], tol=TOL) if i > 0 else 0
        out.append(ker - im)
    return out
