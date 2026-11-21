import os, sys, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import spin_network_ladder, character_model

cache = {}
def ladder(tjmax):
    if tjmax not in cache:
        states, E2, WL, WR = spin_network_ladder(tjmax, lambda *x: 1.0)
        cache[tjmax] = (states, sp.csr_matrix(E2), sp.csr_matrix(WL), sp.csr_matrix(WR))
    return cache[tjmax]

def low(tjmax, g, k=3):
    states, E2, WL, WR = ladder(tjmax)
    n = len(states)
    H = (g**2 / 2) * E2 + (1 / g**2) * (4 * sp.identity(n) - WL - WR)
    vals, vecs = eigsh(H, k=k, which="SA")
    order = np.argsort(vals); vals, vecs = vals[order], vecs[:, order]
    v0 = vecs[:, 0]
    return vals, (v0 @ (E2 @ v0)) / 7, (v0 @ ((WL + WR) @ v0)) / 4      # per link, <tr U_p>/2

# strong-coupling second order (degenerate perturbation theory on the lowest excited doublet)
states, E2, WL, WR = spin_network_ladder(8, lambda *x: 1.0)
E0 = np.diag(E2) / 2; V = -(WL + WR)                    # H = g^2 E0 + (4 + V)/g^2 (dropping constants)
iv = states.index((0, 0, 0)); ex = [states.index((1, 0, 1)), states.index((0, 1, 1))]
def second_order(P, e):
    Q = [k for k in range(len(states)) if k not in P]
    return np.array([[sum(V[i, k] * V[k, j] / (e - E0[k]) for k in Q) for j in P] for i in P])
c0 = second_order([iv], 0.0)[0, 0]; c1 = np.linalg.eigvalsh(second_order(ex, 0.75 * 2))[0]
print(f"strong coupling: Delta = 1.5 g^2 + {c1 - c0:.6f}/g^6 + ...   (first-order splitting of the doublet: {np.linalg.eigvalsh(V[np.ix_(ex, ex)])})")
print("weak coupling (harmonic, gauge-invariant): Delta -> 2*sqrt(3/2) =", round(2 * np.sqrt(1.5), 5), "; single plaquette 2*sqrt(2) =", round(2 * np.sqrt(2), 5))
print()
print("   g    ladder Delta (cutoff j<=10 / j<=14)   strong       1-plaq Delta   <E^2>/link   <tr U_p>/2")
for g in (0.35, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0):
    a, e2a, wa = low(20, g); b, e2b, wb = low(28, g)
    gap_a, gap_b = a[1] - a[0], b[1] - b[0]
    one = character_model(400, g)
    strong = 1.5 * g**2 + (c1 - c0) / g**6
    print(f"  {g:4.2f}   {gap_a:9.5f} / {gap_b:9.5f}           {strong:10.4f}   {one[1]-one[0]:9.5f}     {e2b:8.4f}    {wb:8.4f}")
    assert abs(gap_a - gap_b) < 1e-6                                  # converged in the spin cutoff
    if g >= 2.0:
        assert abs(gap_b - strong) < 1e-3                             # strong-coupling series
    if g <= 0.5:
        assert abs(gap_b - 2 * np.sqrt(1.5)) < 0.1                    # approaching the harmonic limit
