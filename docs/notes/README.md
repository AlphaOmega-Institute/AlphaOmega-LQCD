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

    julia --project=lattice -t auto lattice/test/runtests.jl          # LatticeYM, 1712 tests
    julia docs/notes/check_spectral_transfer.jl                        # 4260 checks
    python3 hamiltonian/check_all.py                                   # NumPy and SciPy
    gunzip -k lattice/results/glueball_b570_L8T16.txt.gz           # 3000 configurations, 2.6 MB
    julia --project=lattice lattice/examples/glueball_analysis.jl lattice/results/glueball_b570_L8T16.txt --bin 10 --ops 2,3,4,5

The data file was produced with

    julia --project=lattice -t 4 lattice/examples/glueball.jl --beta 5.7 --L 8 --T 16 --therm 300 \
          --confs 3000 --skip 4 --levels 0,4,8,16,24,32 --wlevel 16 --rmax 4 --tmax 6 --seed 57

Julia 1.13.4, four threads. Results do not depend on the thread count.

## Results at β = 5.7

8³ × 16, 3000 configurations (4 sweeps apart), bins of 10, GEVP on smearing levels 4, 8, 16, 24 with
t₀ = 0. Full output: `lattice/results/glueball_b570_L8T16_analysis.txt`.

| quantity | value |
|---|---|
| ⟨P⟩ | 0.549239(56) |
| σa² (Cornell fit, e = π/12 / e free) | 0.1677(7) / 0.1576(24) |
| a·m₀⁺⁺ (GEVP, t = 1) | 0.967(33) |
| m₀⁺⁺/√σ | 2.36(8), below the continuum ≈ 3.5: the known scalar dip near β = 5.7 |
| a·E₁ (GEVP, t = 0 → 1) | 2.43(5), poorly determined beyond t = 0 |

Per smearing level (stout steps, ρ = 0.1):

| level | 0 | 4 | 8 | 16 | 24 | 32 | GEVP |
|---|---|---|---|---|---|---|---|
| f, variational | 0.356(7) | 0.896(7) | 0.927(5) | 0.770(10) | 0.600(12) | 0.457(12) | |
| f, from decay at t = 2 | 0.36(3) | 0.83(3) | 0.86(3) | 0.74(3) | 0.57(3) | 0.43(3) | |
| a·δ bound, (E₁ − E₀)√(1/f − 1) | 1.97(8) | 0.50(3) | 0.41(2) | 0.80(4) | 1.20(6) | 1.60(8) | |
| δ_T = ‖(I − P_ψ)Tψ‖ | 0.125(17) | 0.122(17) | 0.133(16) | 0.152(15) | 0.152(15) | 0.145(16) | 0.100(20) |
| window gap 1 − c(1)/c(0) | 0.811(5) | 0.676(5) | 0.673(5) | 0.711(5) | 0.755(5) | 0.798(5) | 0.648(5) |
| δ_T / [(1 − c)/2√2] | 0.44(6) | 0.51(7) | 0.56(7) | 0.60(6) | 0.57(6) | 0.51(6) | 0.43(9) |

Reading the two defect measures together:
- With H, the bound built from f and E₁ puts the best trial state's defect at ≳ 0.41 against a mass gap
  of 0.97, above the Proposition 2 scale m/2√2 ≈ 0.34. The bound depends on the poorly measured E₁.
- With T, the directly measured residual of the best trial states is 0.43–0.56 of the Proposition 2
  scale. Euclidean time damps excited admixtures; this is a necessary condition only (see Status).

## Status

- The β = 5.7 point reproduces known lattice physics: the string tension and the low 0⁺⁺ mass of the
  "scalar dip" near this coupling. It validates the code and the analysis, not a new claim.
- The question the measurements are meant to answer: does the best presentation's defect shrink
  relative to the gap as the lattice gets finer? Corollary 3 of the note needs that, in Wilsonian form
  δ/Δ ~ (am)^p with p > 0. Runs at β = 5.85 and 6.0 are planned. Optimized single-operator overlaps are
  usually close to constant in β, so a flat trend is the expected outcome. That would say single smeared
  operators are too crude a family of presentations, not that the criterion fails.
- The two ground-state-fraction estimators (variational, and from the correlator decay) differ
  systematically. Neither is a bound, and their spread is not an uncertainty interval.
- Which gap goes where. In Proposition 1(ii), η = δ/(Λ − λ): for a one-state window ψ ≈ glueball,
  Λ − λ ≈ E₁ − E₀, the spacing to the next scalar level. In Proposition 2 the tolerance κ < Δ_c compares
  the defect with the coarse gap, which for the window {vacuum, ψ} is the vacuum-to-glueball gap. Since
  ε ≥ δ, κ ≥ 2√2·δ.
- The transfer-matrix residual δ_T = ‖(I − P_ψ)Tψ‖ with T = exp(−aH), from c(2)/c(0) − (c(1)/c(0))²,
  is the defect of the window {vacuum, ψ} with respect to T. It needs no fitted E₁. Proposition 2
  applied to −T compares it with the window gap 1 − c(1)/c(0). This comparison is necessary but not
  sufficient: Proposition 2 also needs the top of T on the window's complement (through η), which three
  correlator values do not determine. A poorly overlapping operator can show a small ratio while the
  glueball itself lies mostly in the complement, and then Proposition 2 does not apply. It is a
  one-lattice residual, not a test of a map between lattice spacings.
- Euclidean time damps excited states, so a defect measured with T is much smaller than the same
  admixture measured with H. A transfer-matrix form of Corollary 3 (maps intertwining T_k² with
  T_{k+1}) is the natural next statement to write down.

## Points for review, most important first

1. Propositions 1 and 2 of the note: the proofs, and the windowed form with W = W₀Q.
2. Corollary 3 and Remark 4: whether the reduction of the gap to summable windowed defects and omitted
   energies at the cutoff is stated correctly, and whether anything essential is missing.
3. Remark 5 and `lattice/examples/glueball_analysis.jl`: the identification of the ground-state fraction
   of a smeared operator with ‖Pφ₀‖² for a one-dimensional window, the defect bound built from it, and
   the transfer-matrix residual (`transfer_residual`) with its use in Proposition 2 for −T.
4. `lattice/src/smearing.jl` and `lattice/src/spectroscopy.jl`: gauge covariance of the stout step, the
   correlator estimators, and the GEVP normalization.
5. `hamiltonian/`: the 6j matrix elements of the spin-network ladder (validated against brute force) and
   the BRST-complex signature (+5, −5).
