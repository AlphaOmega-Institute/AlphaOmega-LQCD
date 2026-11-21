# Small SU(2) Hamiltonian lattice models

Exact computations for the Kogut–Susskind Hamiltonian of SU(2) lattice gauge theory on one plaquette
and on two plaquettes sharing a link:

    H = (g²/2) Σ_links E_l² + (1/g²) Σ_plaquettes (2 − tr U_p).

These models are a first test of the polar formalism on gauge theory: the Gauss law as a BRST
complex, the electric term as d*d, and the gap as the closed-range constant of that map. They are
far too small to say anything about the continuum limit; see "What this does and does not show".

Python 3 with NumPy and SciPy. Run every check with

    python3 hamiltonian/check_all.py

## Files

| File | Contents |
|---|---|
| `su2.py` | Clebsch–Gordan and 6j symbols (Racah formulas, exact rationals), spin matrices, and L²(SU(2)) truncated at spin j_max in the Peter–Weyl basis, with left and right generators and multiplication by the fundamental matrix elements. |
| `models.py` | The single plaquette with all four links and the Gauss law at every vertex; the reduced character model; the two-plaquette ladder by brute force (three link variables, Gauss law at both trivalent nodes) and in a spin-network basis with 6j matrix elements. |
| `brst_chain.py` | The charge of a chain of maps and its cohomology, from the BRST working notes in AlphaOmega-Theory (`manuscripts/notes/brst/test_brst_families.py`). |
| `check_su2.py` | CG and 6j orthogonality and known values; generator algebra; the character identity χ½ χ_j = χ_{j−½} + χ_{j+½}. |
| `plaquette.py` | One link (closed range from compactness); the single plaquette two ways; the Gauss law as a Chevalley–Eilenberg complex with the chain charge; the gap across couplings. |
| `validate_ladder.py`, `validate_ladder2.py` | The spin-network ladder against brute force at j_max = ½, 1 and 3/2. |
| `ladder_gap.py` | The ladder gap across couplings, converged in the spin cutoff, with the strong- and weak-coupling limits and two vacuum observables. |

## Results

**One link.** The electric Laplacian is d*d with d the left-invariant gradient. On SU(2) its spectrum is
j(j+1), so the closed-range constant is c² = 3/4. On a circle it is 1. On a segment of length L it is
about π²/L², which tends to zero: the line has no gap.

**One plaquette.** All four links with an explicit Gauss law (625 states at j_max = ½, 38 416 at
j_max = 1) reproduce the spectrum of the reduced character model exactly. As a Chevalley–Eilenberg
complex for the residual SU(2), the Gauss law has cohomology in ghost numbers 0 and 3, and the chain
charge gives signature (+5, −5) at j_max = 2: the physical states and a negative-norm copy, the familiar
ghost-number issue of BRST for compact groups.

**Two plaquettes.** The spin-network basis agrees with brute force to 10⁻¹³ at j_max = 3/2
(27 000 states, 23 gauge invariant).

**Gap across couplings** (lattice units; converged in the spin cutoff):

| g | two plaquettes | one plaquette | ⟨tr U_p⟩/2 in the vacuum |
|---|---|---|---|
| 0.35 | 2.411 | 2.750 | 0.94 |
| 0.7 | 2.272 | 2.487 | 0.74 |
| 1.0 | 2.105 | 2.196 | 0.47 |
| 1.5 | 3.440 | 3.456 | 0.13 |
| 3.0 | 13.501 | 13.501 | 0.01 |

Strong coupling: Δ = 1.5 g² + 0.933/g⁶ (one plaquette, 14/15) and 1.5 g² + 0.743/g⁶ (two plaquettes).
Weak coupling: the harmonic limits 2√2 and 2√(3/2). The gap is smallest at the crossover from the
electric to the magnetic vacuum.

## What this does and does not show

At strong coupling the gap comes from compactness of SU(2): the electric map has closed range, and four
links times 3/4 times g²/2 give 1.5 g². As the vacuum turns magnetic that protection carries a factor g²
that goes to zero. The weak-coupling gap on these lattices is a finite-volume oscillator gap, and it
already decreases with volume (2.83 for one plaquette, 2.45 for two). Whether a gap survives the
continuum limit is beyond any lattice this small; the glueball measurements in `lattice/` address it
numerically.
