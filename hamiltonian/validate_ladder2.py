import os, sys, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.linalg import orth
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import full_ladder, spin_network_ladder, spin_network_states
tjmax = 3
gauss, E2, WL, WR = full_ladder(tjmax)
G2 = sum(G.conj().T @ G for G in gauss)
nphys = len(spin_network_states(tjmax))
vals, vecs = eigsh(G2, k=nphys + 6, sigma=-0.5, which="LM")
phys = orth(vecs[:, vals < 1e-9])
print(f"jmax={tjmax/2}: brute-force dim {E2.shape[0]}, gauge-invariant {phys.shape[1]}, spin networks {nphys}, next G^2 eigenvalue {np.sort(vals)[nphys]:.3f}")
states, E2s, WLs, WRs = spin_network_ladder(tjmax, lambda *x: 1.0)
for g in (0.6, 1.0, 1.7, 3.0):
    Hb = phys.conj().T @ ((g**2 / 2) * E2 + (1 / g**2) * (4 * sp.identity(E2.shape[0]) - WL - WR)) @ phys
    Hs = (g**2 / 2) * E2s + (1 / g**2) * (4 * np.eye(len(states)) - WLs - WRs)
    eb, es = np.linalg.eigvalsh((Hb + Hb.conj().T) / 2), np.linalg.eigvalsh(Hs)
    print(f"   g={g}: max spectral difference {np.abs(eb - es).max():.1e}")
    assert np.abs(eb - es).max() < 1e-10
assert phys.shape[1] == nphys
