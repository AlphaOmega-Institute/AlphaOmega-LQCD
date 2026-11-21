import os, sys, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.linalg import orth
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from models import full_ladder, spin_network_ladder, spin_network_states
rules = {"+1": lambda a, ap, c, cp, b: 1.0,
         "(-1)^(a+a'+c+c')/2": lambda a, ap, c, cp, b: (-1.0) ** ((a + ap + c + cp) // 2),
         "(-1)^(a'+c'+b+1)/2": lambda a, ap, c, cp, b: (-1.0) ** ((ap + cp + b + 1) // 2),
         "(-1)^(a+c'+b+1)/2": lambda a, ap, c, cp, b: (-1.0) ** ((a + cp + b + 1) // 2)}
for tjmax in (1, 2):
    gauss, E2, WL, WR = full_ladder(tjmax)
    G2 = sum(G.conj().T @ G for G in gauss)
    nphys = len(spin_network_states(tjmax))
    vals, vecs = eigsh(G2, k=nphys + 4, sigma=-0.5, which="LM")
    phys = orth(vecs[:, vals < 1e-9])
    print(f"jmax={tjmax/2}: brute-force dim {E2.shape[0]}, gauge-invariant states {phys.shape[1]}, spin networks {nphys}")
    for name, rule in rules.items():
        states, E2s, WLs, WRs = spin_network_ladder(tjmax, rule)
        ok = True
        for g in (0.7, 1.2, 2.5):
            Hb = phys.conj().T @ ((g**2 / 2) * E2 + (1 / g**2) * (4 * sp.identity(E2.shape[0]) - WL - WR)) @ phys
            Hs = (g**2 / 2) * E2s + (1 / g**2) * (4 * np.eye(len(states)) - WLs - WRs)
            ok &= np.allclose(np.linalg.eigvalsh((Hb + Hb.conj().T) / 2), np.linalg.eigvalsh(Hs), atol=1e-9)
        print(f"   sign rule {name:22s}: spectra agree = {ok}")
        if name == "+1":
            assert ok, "the sign rule used by spin_network_ladder must reproduce the brute-force spectrum"
