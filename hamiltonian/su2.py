"""SU(2) tools for small Hamiltonian lattice gauge models.

Spins are passed doubled (tj = 2j) so that everything is an integer. Clebsch-Gordan and 6j
symbols use the Racah formulas in exact rational arithmetic, then a square root.
"""
from fractions import Fraction
from functools import lru_cache
from math import factorial, sqrt
import numpy as np


def _f(x2):                       # factorial of a half-integer-free quantity given doubled
    assert x2 % 2 == 0 and x2 >= 0, x2
    return factorial(x2 // 2)


@lru_cache(maxsize=None)
def cg(tj1, tm1, tj2, tm2, tJ, tM):
    """<j1 m1 j2 m2 | J M> (Condon-Shortley), doubled arguments."""
    if tm1 + tm2 != tM or not (abs(tj1 - tj2) <= tJ <= tj1 + tj2) or (tj1 + tj2 + tJ) % 2:
        return 0.0
    if abs(tm1) > tj1 or abs(tm2) > tj2 or abs(tM) > tJ:
        return 0.0
    if (tj1 + tm1) % 2 or (tj2 + tm2) % 2 or (tJ + tM) % 2:
        return 0.0
    pre = Fraction((tJ + 1) * _f(tJ + tj1 - tj2) * _f(tJ - tj1 + tj2) * _f(tj1 + tj2 - tJ),
                   _f(tj1 + tj2 + tJ + 2))
    pre *= _f(tJ + tM) * _f(tJ - tM) * _f(tj1 - tm1) * _f(tj1 + tm1) * _f(tj2 - tm2) * _f(tj2 + tm2)
    s = Fraction(0)
    for k in range(0, (tj1 + tj2 + tJ) // 2 + 1):
        args = [2 * k, tj1 + tj2 - tJ - 2 * k, tj1 - tm1 - 2 * k, tj2 + tm2 - 2 * k,
                tJ - tj2 + tm1 + 2 * k, tJ - tj1 - tm2 + 2 * k]
        if min(args) < 0:
            continue
        d = 1
        for a in args:
            d *= _f(a)
        s += Fraction((-1) ** k, d)
    val = float(s) * sqrt(float(pre))
    return val


def _tri(a, b, c):
    return Fraction(_f(a + b - c) * _f(a - b + c) * _f(-a + b + c), _f(a + b + c + 2))


@lru_cache(maxsize=None)
def sixj(a, b, c, d, e, f):
    """{a b c; d e f}, doubled arguments (Racah formula)."""
    for x, y, z in ((a, b, c), (a, e, f), (d, b, f), (d, e, c)):
        if not (abs(x - y) <= z <= x + y) or (x + y + z) % 2:
            return 0.0
    pre = _tri(a, b, c) * _tri(a, e, f) * _tri(d, b, f) * _tri(d, e, c)
    s = Fraction(0)
    lo = max(a + b + c, a + e + f, d + b + f, d + e + c)
    hi = min(a + b + d + e, a + c + d + f, b + c + e + f)
    for t2 in range(lo, hi + 1, 2):
        den = (_f(t2 - a - b - c) * _f(t2 - a - e - f) * _f(t2 - d - b - f) * _f(t2 - d - e - c)
               * _f(a + b + d + e - t2) * _f(a + c + d + f - t2) * _f(b + c + e + f - t2))
        s += Fraction((-1) ** (t2 // 2) * _f(t2 + 2), den)
    return float(s) * sqrt(float(pre))


def generators(tj):
    """Spin-j matrices (Jx, Jy, Jz) in the basis m = j, j-1, ..., -j."""
    j = tj / 2
    ms = [j - k for k in range(tj + 1)]
    Jz = np.diag(ms).astype(complex)
    Jp = np.zeros((tj + 1, tj + 1), complex)
    for k in range(1, tj + 1):
        m = ms[k]
        Jp[k - 1, k] = sqrt(j * (j + 1) - m * (m + 1))
    Jm = Jp.conj().T
    return (Jp + Jm) / 2, (Jp - Jm) / 2j, Jz


class PeterWeyl:
    """L^2(SU(2)) truncated at spin jmax, basis f^j_{mn} = sqrt(2j+1) D^j_{mn}, orthonormal.

    Left regular action (lambda_g f)(U) = f(g^-1 U) has generators L^a = conj(T^a) on the row
    index; right regular action (rho_g f)(U) = f(U g) has generators R^a = T^a on the column
    index. Both square to the Casimir j(j+1).
    """

    def __init__(self, tjmax):
        self.tjmax = tjmax
        self.index = []                                   # (tj, im, in)
        self.blocks = {}
        for tj in range(tjmax + 1):
            start = len(self.index)
            for im in range(tj + 1):
                for in_ in range(tj + 1):
                    self.index.append((tj, im, in_))
            self.blocks[tj] = (start, len(self.index))
        self.dim = len(self.index)

    def left(self):
        out = [np.zeros((self.dim, self.dim), complex) for _ in range(3)]
        for tj, (s, e) in self.blocks.items():
            for a, T in enumerate(generators(tj)):
                out[a][s:e, s:e] = np.kron(T.conj(), np.eye(tj + 1))
        return out

    def right(self):
        out = [np.zeros((self.dim, self.dim), complex) for _ in range(3)]
        for tj, (s, e) in self.blocks.items():
            for a, T in enumerate(generators(tj)):
                out[a][s:e, s:e] = np.kron(np.eye(tj + 1), T)
        return out

    def casimir(self):
        return np.diag([tj / 2 * (tj / 2 + 1) for tj, _, _ in self.index]).astype(complex)

    def fundamental(self, alpha, beta):
        """Multiplication by D^{1/2}_{alpha beta}(U), truncated at jmax (alpha, beta in 0, 1)."""
        ta, tb = 1 - 2 * alpha, 1 - 2 * beta                  # doubled m of the spin-1/2 entries
        M = np.zeros((self.dim, self.dim), complex)
        pos = {k: i for i, k in enumerate(self.index)}
        for i, (tj, im, in_) in enumerate(self.index):
            tm, tn = tj - 2 * im, tj - 2 * in_
            for tJ in (tj - 1, tj + 1):
                if tJ < 0 or tJ > self.tjmax:
                    continue
                tM, tN = tm + ta, tn + tb
                if abs(tM) > tJ or abs(tN) > tJ:
                    continue
                c = cg(1, ta, tj, tm, tJ, tM) * cg(1, tb, tj, tn, tJ, tN)
                if c == 0:
                    continue
                k = pos[(tJ, (tJ - tM) // 2, (tJ - tN) // 2)]
                M[k, i] += sqrt((tj + 1) / (tJ + 1)) * c
        return M

    def trace_fundamental(self):
        return self.fundamental(0, 0) + self.fundamental(1, 1)
