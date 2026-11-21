"""Model builders for step 1."""
import itertools, numpy as np
import scipy.sparse as sp
from su2 import PeterWeyl


def full_plaquette(tjmax):
    pw = PeterWeyl(tjmax); n = pw.dim; I = sp.identity(n, format="csr")
    def on(link, op):
        mats = [I] * 4; mats[link] = sp.csr_matrix(op)
        out = mats[0]
        for m in mats[1:]:
            out = sp.kron(out, m, format="csr")
        return out
    L, R, C = pw.left(), pw.right(), pw.casimir()
    # links l0..l3 oriented around the plaquette, U_p = U0 U1 U2 U3; vertex v joins end of l(v-1) and start of l(v)
    gauss = [on((v - 1) % 4, R[a]) - on(v, L[a]) for v in range(4) for a in range(3)]
    E2 = sum(on(l, C) for l in range(4))
    F = {(a, b): pw.fundamental(a, b) for a in (0, 1) for b in (0, 1)}
    W = None
    for a, b, c, e in itertools.product((0, 1), repeat=4):
        term = sp.kron(sp.kron(sp.kron(sp.csr_matrix(F[a, b]), sp.csr_matrix(F[b, c])), sp.csr_matrix(F[c, e])), sp.csr_matrix(F[e, a]), format="csr")
        W = term if W is None else W + term
    return pw, gauss, E2, W

def character_model(tjmax, g):
    tj = np.arange(tjmax + 1); j = tj / 2
    H = np.diag(2 * g**2 * j * (j + 1) + 2 / g**2)
    H -= (np.eye(tjmax + 1, k=1) + np.eye(tjmax + 1, k=-1)) / g**2
    return np.linalg.eigvalsh(H)



# ---- two plaquettes sharing a link (open ladder) ------------------------------------------------
# After the 2-valent vertices are gauge-fixed the ladder is a theta graph: three edges from node 1
# to node 4 -- the left path (3 links, spin a), the right path (3 links, spin b) and the shared rung
# (spin c). Plaquettes: W_L = tr(U_a U_c^-1), W_R = tr(U_b U_c^-1).
LINK_WEIGHTS = (3, 3, 1)                                  # links carried by edges a, b, c


def full_ladder(tjmax):
    """Brute force: three link variables, explicit Gauss law at both trivalent nodes."""
    pw = PeterWeyl(tjmax); n = pw.dim; I = sp.identity(n, format="csr")
    def on(k, op):
        mats = [I] * 3; mats[k] = sp.csr_matrix(op)
        return sp.kron(sp.kron(mats[0], mats[1], format="csr"), mats[2], format="csr")
    L, R, C = pw.left(), pw.right(), pw.casimir()
    gauss = [on(0, L[a]) + on(1, L[a]) + on(2, L[a]) for a in range(3)] + \
            [on(0, R[a]) + on(1, R[a]) + on(2, R[a]) for a in range(3)]
    E2 = sum(w * on(k, C) for k, w in enumerate(LINK_WEIGHTS))
    F = {(a, b): sp.csr_matrix(pw.fundamental(a, b)) for a in (0, 1) for b in (0, 1)}
    WL = sum(sp.kron(sp.kron(F[a, b], I), F[a, b].conj().T, format="csr") for a in (0, 1) for b in (0, 1))
    WR = sum(sp.kron(sp.kron(I, F[a, b]), F[a, b].conj().T, format="csr") for a in (0, 1) for b in (0, 1))
    return gauss, E2, WL, WR


def spin_network_states(tjmax):
    return [(a, b, c) for a in range(tjmax + 1) for b in range(tjmax + 1) for c in range(tjmax + 1)
            if abs(a - b) <= c <= a + b and (a + b + c) % 2 == 0]


def spin_network_ladder(tjmax, sign_rule):
    """Spin-network basis |a b c> (doubled spins). Electric term diagonal; plaquettes by 6j symbols."""
    from su2 import sixj
    states = spin_network_states(tjmax)
    pos = {s: i for i, s in enumerate(states)}
    n = len(states)
    E2 = np.diag([sum(w * (t / 2) * (t / 2 + 1) for w, t in zip(LINK_WEIGHTS, s)) for s in states])
    WL = np.zeros((n, n)); WR = np.zeros((n, n))
    for (a, b, c), i in pos.items():
        for da in (-1, 1):
            for dc in (-1, 1):
                ap, cp = a + da, c + dc
                if (ap, b, cp) in pos:          # left plaquette moves a and c
                    val = np.sqrt((a + 1) * (ap + 1) * (c + 1) * (cp + 1)) * sixj(a, c, b, cp, ap, 1) ** 2
                    WL[pos[(ap, b, cp)], i] = sign_rule(a, ap, c, cp, b) * val
                bp = b + da
                if (a, bp, cp) in pos:          # right plaquette moves b and c
                    val = np.sqrt((b + 1) * (bp + 1) * (c + 1) * (cp + 1)) * sixj(b, c, a, cp, bp, 1) ** 2
                    WR[pos[(a, bp, cp)], i] = sign_rule(b, bp, c, cp, a) * val
    return states, E2, WL, WR
