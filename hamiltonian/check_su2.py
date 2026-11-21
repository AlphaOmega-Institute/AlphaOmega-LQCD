import os, sys, numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from su2 import cg, sixj, generators, PeterWeyl
# CG: orthogonality and known values
assert abs(cg(1, 1, 1, -1, 0, 0) - 1/np.sqrt(2)) < 1e-14 and abs(cg(1, -1, 1, 1, 0, 0) + 1/np.sqrt(2)) < 1e-14
for tj1, tj2 in ((1, 1), (2, 1), (3, 2), (4, 3)):
    rows = []
    for tJ in range(abs(tj1 - tj2), tj1 + tj2 + 1, 2):
        for tM in range(-tJ, tJ + 1, 2):
            rows.append([cg(tj1, tm1, tj2, tm2, tJ, tM) for tm1 in range(-tj1, tj1 + 1, 2) for tm2 in range(-tj2, tj2 + 1, 2)])
    U = np.array(rows)
    assert np.allclose(U @ U.T, np.eye(len(U)))
# 6j: known values and orthogonality
assert abs(sixj(1, 1, 0, 0, 0, 1) + 1/np.sqrt(2)) < 1e-14                     # {1/2 1/2 0; 0 0 1/2}
for a, b, d, e in ((2, 2, 2, 2), (3, 1, 2, 2), (4, 2, 3, 1)):
    for f in range(0, 9):
        for fp in range(0, 9):
            s = sum((c + 1) * (f + 1) * sixj(a, b, c, d, e, f) * sixj(a, b, c, d, e, fp) for c in range(0, 9))
            ok = all(abs(x - y) <= z <= x + y and (x + y + z) % 2 == 0 for x, y, z in ((a, e, f), (d, b, f)))
            assert abs(s - (1.0 if (f == fp and ok) else 0.0)) < 1e-12, (a, b, d, e, f, fp, s)
# generators
for tj in range(1, 5):
    X, Y, Z = generators(tj)
    assert np.allclose(X @ Y - Y @ X, 1j * Z)
# Peter-Weyl: L and R are commuting representations (L anti-), both square to the Casimir
pw = PeterWeyl(4)
L, R, C = pw.left(), pw.right(), pw.casimir()
assert np.allclose(sum(x @ x for x in L), C) and np.allclose(sum(x @ x for x in R), C)
assert all(np.allclose(l @ r, r @ l) for l in L for r in R)
assert np.allclose(R[0] @ R[1] - R[1] @ R[0], 1j * R[2])
assert np.allclose(-L[0] @ -L[1] - (-L[1]) @ -L[0], 1j * -L[2])
# multiplication by the fundamental: chi_1/2 chi_j = chi_{j-1/2} + chi_{j+1/2} on class functions
chi = pw.trace_fundamental()
def character(tj):
    v = np.zeros(pw.dim, complex)
    s, e = pw.blocks[tj]
    for i in range(s, e):
        _, im, in_ = pw.index[i]
        if im == in_:
            v[i] = 1 / np.sqrt(tj + 1)                 # chi_j = sum_m D_mm = sum_m f_mm / sqrt(2j+1)
    return v
for tj in range(0, 4):
    w = chi @ character(tj)
    expect = (character(tj - 1) if tj >= 1 else 0) + character(tj + 1)
    assert np.allclose(w, expect), tj
# multiplication operators commute with each other (functions) and fundamental entries are covariant
F = {(a, b): pw.fundamental(a, b) for a in (0, 1) for b in (0, 1)}
inner = np.ix_(*(2 * [[i for i, (tj, _, _) in enumerate(pw.index) if tj <= 2]]))
assert np.allclose((F[0, 0] @ F[1, 1])[inner], (F[1, 1] @ F[0, 0])[inner])
print("su2 checks passed")
