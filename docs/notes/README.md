# Working notes

| File | Contents |
|---|---|
| `spectral_transfer.tex` / `.pdf` | *Spectral transfer along block-spin maps.* A block-spin map written as a partial isometry; coarse and fine levels shadow each other within its defect; a coarse gap transfers to the fine lattice; the renormalization-group chain this implies; and (Proposition 6) a bound on the defect of a map between two lattices compared through operator windows. |
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

    julia --project=lattice -t auto lattice/test/runtests.jl          # LatticeYM, 2218 tests
    julia docs/notes/check_spectral_transfer.jl                        # 6563 checks
    python3 hamiltonian/check_all.py                                   # NumPy and SciPy
    gunzip -k lattice/results/glueball_b570_L8T16.txt.gz           # 3000 configurations, 2.6 MB
    julia --project=lattice lattice/examples/glueball_analysis.jl lattice/results/glueball_b570_L8T16.txt \
          --bin 10 --ops 2,3,4,5 --windows "2,3,4,5;2,3"

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
| σa², Cornell fit with e = π/12 fixed | 0.1677(7) |
| σa², Cornell fit with e free (e = 0.342(16)) | 0.1576(24) |
| a·m₀⁺⁺ (GEVP, t = 1) | 0.967(33) |
| m₀⁺⁺/√σ (with the fixed-e σ) | 2.36(8), below the continuum ≈ 3.5: the known scalar dip near β = 5.7 |
| a·E₁ (GEVP, t = 0 → 1) | 2.43(5), poorly determined beyond t = 0 |

Comparison with Lucini and Teper (hep-lat/0103027, Table 2; β = 5.70 on 8⁴, string tension from
flux loops): a√σ = 0.3961(23), i.e. σa² = 0.1569(18); am₀⁺⁺ = 0.999(25); first excited scalar
am₀⁺⁺⋆ = 2.07(21). The free-e fit agrees with their string tension and the fixed-e fit is 7% higher, a
systematic spread larger than either statistical error. Our scalar mass agrees. Our E₁ is probably
contaminated upward.

Per smearing level (stout steps, ρ = 0.1):

| level | 0 | 4 | 8 | 16 | 24 | 32 | GEVP |
|---|---|---|---|---|---|---|---|
| f, variational | 0.356(7) | 0.896(7) | 0.927(5) | 0.770(10) | 0.600(12) | 0.457(12) | |
| f, from decay at t = 2 | 0.36(3) | 0.83(3) | 0.86(3) | 0.74(3) | 0.57(3) | 0.43(3) | |
| a·δ bound, (E₁ − E₀)√(f(1 − f)) | 0.70(3) | 0.45(2) | 0.38(2) | 0.62(3) | 0.72(3) | 0.73(3) | |
| δ_T = ‖(I − P_ψ)Tψ‖ | 0.125(17) | 0.122(17) | 0.133(16) | 0.152(15) | 0.152(15) | 0.145(16) | 0.100(20) |
| window gap 1 − c(1)/c(0) | 0.811(5) | 0.676(5) | 0.673(5) | 0.711(5) | 0.755(5) | 0.798(5) | 0.648(5) |
| δ_T / [(1 − c)/2√2] | 0.44(6) | 0.51(7) | 0.56(7) | 0.60(6) | 0.57(6) | 0.51(6) | 0.43(9) |

Reading the two defect measures together:
- With H, the sharp bound given f and the gap E₁ − E₀ (Remark 5) puts the best trial state's residual
  at ≳ 0.38, or ≳ 0.29 with the published E₁. Both inputs are estimates, so neither figure establishes
  that the gap-transfer criterion fails. (An earlier version substituted E₁ for the bottom of the
  omitted-space compression and overstated the bound by 1/f.)
- With T, the directly measured residual of the best trial states is 0.43–0.56 of the Proposition 2
  scale. Euclidean time damps excited admixtures; this is a necessary condition only (see Status).

## Two lattices: pilot on the coarse ensemble

Proposition 6 compares the window compression K and leakage L of a coarse lattice (one time step) with
those of a blocked fine lattice (two time steps). Before running the fine lattice, the coarse data fix
how precisely these quantities can be measured. Coarse leakage, signed eigenvalues with jackknife errors
(the exact matrix is positive):

| window (stout levels) | eigenvalues of L_c | largest | eigenvalues of K_c |
|---|---|---|---|
| 4, 8, 16, 24 | −0.0100(50), −0.0071(44), 0.0038(44), 0.0111(42) | 0.0111(42) | 0.352(5), 0.088(5), 0.024(5), 0.005(5) |
| 4, 8 | 0.0039(44), 0.0149(42) | 0.0149(42) | 0.339(5), 0.064(5) |

Noise floor of the in-window mismatch (`examples/window_noise_floor.jl`): pairs of independent
block-bootstrap resamples of 3000 configurations from this ensemble, whose true mismatch is zero, give an
apparent mismatch whose 95th percentile is 0.017 for levels 4, 8 and 0.022 for levels 4, 8, 16, 24 (2000
pairs; the review found the same with 4000). Levels 4 and 8 are the primary window, all four a
sensitivity check.

The negative eigenvalues of the larger window are noise of the same size as its positive ones. The
largest eigenvalue of a noisy matrix is also biased upward, so its square root, 0.105(20) or 0.122(17),
is an estimate dominated by noise, not a bound. Two standard errors above it gives 0.14–0.15.

Budget. In transfer-matrix units the gap of the window {vacuum, glueball} is 1 − 0.35 ≈ 0.65, so the
Proposition 2 scale is about 0.65/2√2 ≈ 0.23. With leakage up to 0.15, an in-window mismatch below about
0.08 keeps even the upper bound of Proposition 6 under that scale; this budget is conservative, because
the analysis also computes the defect itself (the equality in Proposition 6), which can stay under the
scale with a larger mismatch. The
eigenvalues of K_c are measured to ±0.005, and the fine lattice, at the same physical volume and with
lags 0, 2, 4, should reach similar precision with a similar number of configurations. So the fine run can
resolve the comparison if the mismatch is above about 0.02.

Fine run (to be done on a faster machine): β = 6.054, where the Necco–Sommer interpolation gives
a(5.7)/a(6.054) = 2.0015; 16³ × 32 blocked once to 8³ × 32; blocking parameters α = 0, 0.3, 0.6 on the
same configurations; stout levels 4, 8, 16, 24. Each α is measured in two variants: the fixed blocking
origin, and the zero-momentum projection (operator amplitudes averaged over the eight origins before the
correlators are formed). Target 3000 configurations, with the first 1000 inspected as a pilot while the
same chain continues; `--checkpoint`/`--resume` make an interrupted run continue exactly. Cost:
the measurement now smears 24 blocked lattices per configuration (3 α × 8 origins) up to 24 stout
steps, which costs more than the sweeps: well over 13 s per configuration on the four threads used here
(being re-timed). Time it first on the machine that runs it, with `--confs 25`.

## Status

- The β = 5.7 point reproduces known lattice physics: the string tension and the low 0⁺⁺ mass of the
  "scalar dip" near this coupling. It validates the code and the analysis, not a new claim.
- Corrections after review: Remark 5 now uses the bottom of the omitted-space compression, with the sharp
  bound δ ≥ g√(f(1 − f)); Corollary 3 requires each window to contain the simple coarse vacuum and the
  first excited level, and the continuum statement a uniform margin; the variational overlap estimator
  converts its amplitude to t = 0 with exp(E₀t₀) (it used exp(2E₀t₀); the committed run used t₀ = 0 and
  is unaffected). Blocking uses a polar SU(3) projection, covariant at both ends of a path.
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

New in this round:

1. Proposition 6 and Remark 7: the bound and the equality for the defect itself, the alignment induced
   by blocking as the polar part of the operator-state correspondence, centered excitations, and the
   fixed-origin and zero-momentum variants of the blocked operators.
2. `block_spatial` (now with an origin) and `su3_project` (with a covariant fallback, the straight path,
   for numerically singular path sums), `window_defect`, and the dimension and row-width checks of
   `window_transfer.jl`.
3. The pilot and budget above: whether the fine run, as designed, can answer the question.
4. `examples/blocked_glueball.jl` and `examples/window_transfer.jl`: the two-ensemble jackknife and the
   use of fine lag 2 against coarse lag 1.

From the previous round:

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
