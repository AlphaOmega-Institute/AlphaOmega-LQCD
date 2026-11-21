# Working notes

| File | Contents |
|---|---|
| `spectral_transfer.tex` / `.pdf` | *Spectral transfer along block-spin maps.* A block-spin map written as a partial isometry; coarse and fine levels shadow each other within its defect; a coarse gap transfers to the fine lattice; the renormalization-group chain this implies. |
| `check_spectral_transfer.jl` | Random-instance checks of every inequality in the note (Julia, standard library only). |

## Context

These notes test whether the polar formalism of the A|Ω⟩ program gives any leverage on the Yang–Mills
mass gap. Three pieces of work belong together:

1. `hamiltonian/`: exact SU(2) Hamiltonian models on one and two plaquettes. The Gauss law as a BRST
   complex, the gap as a closed-range constant, the strong- and weak-coupling limits.
2. `docs/notes/spectral_transfer.tex`: the transfer criterion. It reduces a gap in physical units to
   block-spin maps with summable low-energy defects and omitted states at the cutoff. It does not
   construct those maps; that is the renormalization-group problem itself.
3. `lattice/`: glueball spectroscopy in four-dimensional SU(3), measuring the lattice counterparts of
   the note's quantities. The state a smeared operator creates is a one-dimensional presentation of the
   glueball, and its ground-state fraction f bounds the defect: δ ≳ (E₁ − E₀)·√(1/f − 1).

## Reproduce

    julia --project=lattice -t auto lattice/test/runtests.jl          # LatticeYM, 1692 tests
    julia docs/notes/check_spectral_transfer.jl                        # 4260 checks
    python3 hamiltonian/check_all.py                                   # NumPy and SciPy
    gunzip -k lattice/data/glueball_b570_L8T16.txt.gz
    julia --project=lattice lattice/examples/glueball_analysis.jl lattice/data/glueball_b570_L8T16.txt --bin 10 --ops 2,3,4,5

The data file was produced with

    julia --project=lattice -t 4 lattice/examples/glueball.jl --beta 5.7 --L 8 --T 16 --therm 300 \
          --confs 3000 --skip 4 --levels 0,4,8,16,24,32 --wlevel 16 --rmax 4 --tmax 6 --seed 57

Julia 1.13.4, four threads. Results do not depend on the thread count.

## Results at β = 5.7

The 3000-configuration run is in progress. Its data file and these results will be added in the next commit; until then the data file named above does not exist yet.

## Status

- The β = 5.7 point reproduces known lattice physics: the string tension and the low 0⁺⁺ mass of the
  "scalar dip" near this coupling. It validates the code and the analysis, not a new claim.
- The question the measurements are meant to answer: does the best presentation's defect shrink
  relative to the gap as the lattice gets finer? Corollary 3 of the note needs that, in Wilsonian form
  δ/Δ ~ (am)^p with p > 0. Runs at β = 5.85 and 6.0 are planned. Optimized single-operator overlaps are
  usually close to constant in β, so a flat trend is the expected outcome. That would say single smeared
  operators are too crude a family of presentations, not that the criterion fails.
- The two ground-state-fraction estimators disagree by several standard errors (see Results). The
  defect bounds are quoted as ranges for that reason.

## Points for review, most important first

1. Propositions 1 and 2 of the note: the proofs, and the windowed form with W = W₀Q.
2. Corollary 3 and Remark 4: whether the reduction of the gap to summable windowed defects and omitted
   energies at the cutoff is stated correctly, and whether anything essential is missing.
3. Remark 5 and `lattice/examples/glueball_analysis.jl`: the identification of the ground-state fraction
   of a smeared operator with ‖Pφ₀‖² for a one-dimensional window, and the defect bound built from it.
4. `lattice/src/smearing.jl` and `lattice/src/spectroscopy.jl`: gauge covariance of the stout step, the
   correlator estimators, and the GEVP normalization.
5. `hamiltonian/`: the 6j matrix elements of the spin-network ladder (validated against brute force) and
   the BRST-complex signature (+5, −5).
