"""Step 1a: single link and single plaquette, SU(2), Kogut-Susskind Hamiltonian
H = (g^2/2) sum_links E_l^2 + (1/g^2) sum_plaquettes (2 - tr U_p)."""
import os, sys, itertools, numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.linalg import orth
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from su2 import PeterWeyl
from brst_chain import chain_charge, cohomology, chain_cohomology_dims
np.set_printoptions(precision=5, suppress=True, linewidth=120)

# ---- 1. one link: the electric map d (left-invariant gradient) and its closed range -------------
print("1. One link: electric Laplacian = d*d, closed-range constant c^2 = smallest nonzero eigenvalue")
pw = PeterWeyl(8)
d = np.vstack(pw.left())                                  # d: L^2 -> L^2 (x) C^3
lap = d.conj().T @ d
w = np.linalg.eigvalsh(lap)
print("   SU(2): spectrum j(j+1) =", np.unique(np.round(w, 10))[:5], "-> c^2 =", np.unique(np.round(w, 10))[1])
print("   U(1) (circle): -d^2/dtheta^2 has spectrum n^2 -> c^2 = 1")
for L in (10.0, 100.0, 1000.0):
    n = 400; x = np.linspace(0, L, n, endpoint=False); h = L / n           # line segment, Dirichlet
    D = (np.eye(n) - np.eye(n, k=-1))[:, :-1] / h                          # gradient
    print(f"   segment of length {L:6.0f}: c^2 = {np.linalg.eigvalsh(D.T @ D)[0]:.2e} (-> 0: the line has no gap)")

# ---- 2. single plaquette with all four links and an explicit Gauss law ---------------------------
from models import full_plaquette, character_model

print("\n2. Single plaquette, four links, explicit Gauss law at the four vertices")
for tjmax in (1, 2):
    pw, gauss, E2, W = full_plaquette(tjmax)
    comm = max(abs((W @ G - G @ W)).max() for G in gauss)
    G2 = sum(G.conj().T @ G for G in gauss)
    vals, vecs = eigsh(G2, k=6, sigma=-0.5, which="LM")
    phys = orth(vecs[:, vals < 1e-9])                       # eigsh vectors need not be orthonormal
    print(f"   jmax={tjmax/2}: dim {E2.shape[0]}, [W, Gauss] max = {comm:.1e}, gauge-invariant states: {phys.shape[1]} (expected {tjmax + 1})")
    for g in (0.8, 1.5, 3.0):
        Hp = phys.conj().T @ ((g**2 / 2) * E2 + (1 / g**2) * (2 * sp.identity(E2.shape[0]) - W)) @ phys
        e_full = np.linalg.eigvalsh((Hp + Hp.conj().T) / 2)
        assert np.allclose(e_full, character_model(tjmax, g), atol=1e-9)
    print(f"      spectrum on the physical states equals the character model at the same truncation")

# ---- 3. Gauss law as BRST: Chevalley-Eilenberg complex for the residual conjugation, our chain charge
print("\n3. Gauss law as a BRST complex (residual SU(2) at the base point), chain charge of the families note")
def ce_complex(X):
    """C^k = V (x) Lambda^k(su2*), d = sum_a rho(e_a) (x) c^a + structure-constant term, rho(e_a) = i X_a."""
    dimV = X[0].shape[0]
    forms = {k: list(itertools.combinations(range(3), k)) for k in range(4)}
    eps = np.zeros((3, 3, 3))
    for a, b, c in itertools.permutations(range(3)):
        eps[a, b, c] = np.linalg.det(np.eye(3)[[a, b, c]])
    f = -eps                                            # [e_b, e_c] = f[a,b,c] e_a with e = iX
    def wedge(a, I):                                    # c^a wedge c^I -> (sign, J)
        if a in I:
            return 0, None
        J = tuple(sorted(I + (a,)))
        return (-1) ** sum(1 for i in I if i < a), J
    maps = []
    for k in range(3):
        src, dst = forms[k], forms[k + 1]
        M = np.zeros((dimV * len(dst), dimV * len(src)), complex)
        for si, I in enumerate(src):
            for a in range(3):                          # rho(e_a) (x) c^a wedge
                s, J = wedge(a, I)
                if s:
                    di = dst.index(J)
                    M[di * dimV:(di + 1) * dimV, si * dimV:(si + 1) * dimV] += s * 1j * X[a]
            for r, a in enumerate(I):                   # derivation: d c^a = -1/2 f[a,b,c] c^b c^c
                rest = I[:r] + I[r + 1:]
                for b, c in itertools.permutations(range(3), 2):
                    coef = -0.5 * f[a, b, c] * (-1) ** r
                    if coef == 0:
                        continue
                    s1, J1 = wedge(c, rest[:r] if False else ())
                    # build c^{I[:r]} c^b c^c c^{I[r+1:]} and sort
                    seq = list(I[:r]) + [b, c] + list(I[r + 1:])
                    if len(set(seq)) < len(seq):
                        continue
                    perm_sign = np.linalg.det(np.eye(len(seq))[np.argsort(seq)]) if seq else 1
                    J = tuple(sorted(seq)); di = dst.index(J)
                    M[di * dimV:(di + 1) * dimV, si * dimV:(si + 1) * dimV] += coef * perm_sign * np.eye(dimV)
        maps.append(M)
    dims = [dimV * len(forms[k]) for k in range(4)]
    return dims, maps
pw = PeterWeyl(4)
X = [r - l for r, l in zip(pw.right(), pw.left())]      # conjugation U -> g U g^-1
dims, maps = ce_complex(X)
assert all(np.allclose(maps[k + 1] @ maps[k], 0) for k in range(2))
print("   d^2 = 0: True")
print("   chain cohomology dims by ghost number 0..3:", chain_cohomology_dims(dims, maps), "(class functions up to j=2: 5)")
Q, J, N = chain_charge(dims, maps, [0, 1, 2, 3])
_, plus, minus = cohomology(Q, J)
print(f"   reduced space of the chain charge: signature (+{plus}, -{minus})")
assert (plus, minus) == (5, 5) and list(chain_cohomology_dims(dims, maps)) == [5, 0, 0, 5]

# ---- 4. gap versus coupling, with strong- and weak-coupling limits ------------------------------
print("\n4. Single-plaquette gap Delta(g) on the physical states (character model, converged in jmax)")
print("     g      Delta      strong 1.5g^2+(14/15)/g^6     weak 2*sqrt(2)")
for g in (0.2, 0.3, 0.5, 0.8, 1.0, 1.3, 2.0, 3.0, 5.0):
    e = character_model(600, g); e2 = character_model(800, g)
    assert abs((e[1] - e[0]) - (e2[1] - e2[0])) < 1e-9
    print(f"   {g:4.1f}  {e[1]-e[0]:9.5f}   {1.5*g**2 + (14/15)/g**6:12.5f}            {2*np.sqrt(2):.5f}")
    if g >= 2.0:
        assert abs((e[1] - e[0]) - (1.5 * g**2 + (14 / 15) / g**6)) < 1e-3
    if g <= 0.3:
        assert abs((e[1] - e[0]) - 2 * np.sqrt(2)) < 0.06
