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
   glueball; if the rest of the channel lies at least g above the glueball, its ground-state fraction f
   bounds the defect from below, sharply: δ ≥ g·√(f(1 − f)) (Remark 5).

## Reproduce

    julia --project=lattice -t auto lattice/test/runtests.jl          # LatticeYM, 3909 tests
    julia --project=lattice lattice/test/runner_checks.jl              # the fine-run program, 19 checks
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
sensitivity check. For the 1000-configuration pilot against the full coarse ensemble (`--size2 1000`)
the 95th percentile for levels 4 and 8 is 0.024, so a pilot mismatch near 0.02 is consistent with
noise alone.

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

Fine run: β = 6.054, where the Necco–Sommer interpolation gives
a(5.7)/a(6.054) = 2.0015; 16³ × 32 blocked once to 8³ × 32; blocking parameters α = 0, 0.3, 0.6 on the
same configurations; stout levels 4, 8, 16, 24. Each α is measured in two variants: the fixed blocking
origin, and the zero-momentum projection (operator amplitudes averaged over the eight origins before the
correlators are formed). Target 3000 configurations, with the first 1000 inspected as a pilot while the
same chain continues. Cost: the
measurement smears 24 blocked lattices per configuration (3 α × 8 origins) up to 24 stout steps, which
costs more than the sweeps: 20 s per configuration on the four threads used here while other jobs ran
(perhaps 13–15 s alone), against 4.5 s for the fixed-origin measurement alone. Time it first on the
machine that runs it, with `--confs 25`. `--levels 4,8` (the primary window only) roughly halves the
cost.

Protecting the measurements (after review of a9dcaad):

- A checkpoint, written every 25 configurations and at the last one, records the length of the output
  written so far and its FNV-1a hash, besides the fields and the random state; the output is flushed
  and synced to disk first. On `--resume` the output must begin with exactly those bytes (same header,
  same rows), so a checkpoint cannot be joined to another run's file or to a modified one; anything
  after them is cut off in place, without rewriting the committed rows, and the run appends. Physics
  options come from the checkpoint, and any given again must agree. A fresh run refuses to overwrite an
  existing output or checkpoint file. `--checkpoint` and `--resume` require `--out`.
- The analysis scripts read only lines ending in a newline (`read_measurements`), check every row read
  (all tokens numbers, constant width, no second header), and ignore an unfinished last line, so the
  file can be analysed while it is written. The writer now emits each row in one write.
- `su3_project` returns the straight path when det of the polar part is within 10⁻⁸ of −1 in phase,
  where the principal cube root jumps by a center element and roundoff picks the side; a configuration
  built to put a path sum on the cut is now blocked covariantly (without the band 17 of 40 random gauge
  transformations changed the blocked link by √3 in norm). In a sample of 983 040 path sums from
  thermalized β = 6.054 configurations (8³ × 16) none reached the band: |arg det| stayed below 1.39 at
  α = 0.3 and below 3.02 at α = 0.6. That does not bound the rate over the 10⁹ path sums of the run;
  the progress line and the checkpoint carry the count of fallbacks, and a nonzero count means the
  safeguard was used.
- `lattice/test/runner_checks.jl` runs the program on small lattices: two hard kills and resumes give
  output byte-identical to an uninterrupted chain; a seed-10 checkpoint is refused on a seed-11 file,
  as are a file with a changed committed row, a truncated file and a conflicting physics option; and
  reads of the growing file always return rows of the finished one.

Running it. In a clone, make a separate run directory at the reviewed commit (later updates to the
branch then do not touch it), time it, run the pilot to a checkpoint boundary, and extend the same
chain if the pilot warrants it:

    git fetch origin <branch>
    git worktree add --detach ../lqcd-fine-run <commit>
    cd ../lqcd-fine-run
    julia --project=lattice -t N lattice/examples/blocked_glueball.jl --confs 25 --therm 2 > /dev/null
    caffeinate -i julia --project=lattice -t N lattice/examples/blocked_glueball.jl --confs 1000 \
          --out fine_b6054.txt --checkpoint fine_b6054.ckpt               # pilot (caffeinate: macOS)
    caffeinate -i julia --project=lattice -t N lattice/examples/blocked_glueball.jl --confs 3000 \
          --resume fine_b6054.ckpt --out fine_b6054.txt                   # after an interruption, or to extend

with N the thread count that timed fastest.

Pilot analysis, at any time during the run:

    gunzip -k lattice/results/glueball_b570_L8T16.txt.gz              # once
    julia --project=lattice lattice/examples/window_transfer.jl lattice/results/glueball_b570_L8T16.txt \
          fine_b6054.txt --window 4,8

and the glueball mass on both lattices in the same physical units (a spectral difference between
the windows bounds the mismatch from below for every alignment, by Weyl's inequality):

    julia --project=lattice lattice/examples/blocked_masses.jl fine_b6054.txt \
          --coarse lattice/results/glueball_b570_L8T16.txt

`window_transfer.jl` ends with the fixed-origin minus zero-momentum differences of the mismatch, the bound, the defect and
the departure from the aligned graph (below) for each α, with paired errors (same configurations and
bins): how much the choice of blocking grid changes the conclusion.

## Pilot at β = 6.054: 1000 configurations

The fine run of 16³ × 32 at β = 6.054 (4 sweeps apart, 12 threads on an Apple M4 Max, 4.8 s per
configuration, no projection fallbacks) reached its first 1000 configurations; the same chain is being
extended to 3000. The output file will be added to `lattice/results/` with the full run. Coarse side:
the committed β = 5.7 ensemble. Window: stout levels 4 and 8; bins of 10; jackknife errors over both
ensembles.

| α, blocking | ‖M − I‖ | ‖K_f V − V K_c‖ | spectral floor | defect | total departure |
|---|---|---|---|---|---|
| 0, fixed origin | 1.735(13) | 0.303(7) | 0.292 | 0.331(11) | 0.180(10) |
| 0, zero momentum | 0.284(8) | 0.253(7) | 0.196 | 0.292(12) | 0.172(10) |
| 0.3, fixed origin | 0.275(7) | 0.210(8) | 0.175 | 0.259(13) | 0.159(12) |
| 0.3, zero momentum | 0.183(5) | 0.174(8) | 0.113 | 0.223(15) | 0.143(13) |
| 0.6, fixed origin | 0.143(6) | 0.195(8) | 0.154 | 0.246(14) | 0.152(12) |
| 0.6, zero momentum | 0.249(4) | 0.167(8) | 0.106 | 0.218(15) | 0.140(14) |
| coarse with itself (control) | 0 | 0 | 0 | 0.122 | 0.122 |

The spectral floor is the largest difference between matched eigenvalues of K_f and K_c; by Weyl's
inequality ‖V*K_f V − K_c‖ cannot fall below it for any unitary V.

- The in-window mismatch is far above the noise reference (0.024) for every recipe, and the window
  {4, 8, 16, 24} gives the same picture (best 0.167(8), floor 0.087). The fine and coarse operators
  that blocking pairs do not evolve alike over one matched interval, and no change of the pairing can
  make them.
- Most of the mismatch is spectral. The top eigenvalues, 0.339(5) coarse and 0.233(7) fine (α = 0.6,
  zero momentum), are close to the decays of the GEVP ground-state correlators over the matched
  interval, but not identical to them: the window maximizes the correlator at fine lag 2, the GEVP of
  the mass analysis at lag 1. With periodic cosh ratios the lag-one combination decays by 0.3393 over
  coarse t = 0 → 1 and by 0.2329 over fine t = 0 → 2 (window: 0.3392, 0.2331). At α = 0, zero
  momentum, it gives 0.1399 against the window's 0.1436: for the better recipe the combination chosen
  at one fine step is nearly optimal at two, for the cruder one it is not.
- Averaging over the eight blocking origins lowers the defect at every α, by 0.038(9), 0.036(9) and
  0.029(8) (paired errors). α = 0.6 with zero momentum is the primary recipe; α = 0.3 is close enough
  that a paired comparison of the two is needed before preferring either.
- Neutrality and invariance rank the recipes in opposite order at α = 0.6: the fixed origin gives the
  more nearly isometric raw map (‖M − I‖ 0.143 against 0.249), while the zero-momentum average tracks
  the evolution better inside the window (mismatch 0.195 against 0.167, paired difference 0.0282(42),
  6.7 errors). In the total departure from the aligned graph, which includes leakage, the improvement
  is 0.0122(69), 1.8 errors.
- The pilot resolves disagreement inside the windows. Out of the windows, the evolution (≈ 0.12) is
  compatible in size with the coarse self-comparison; similar magnitudes leave open which states the
  evolution reaches outside each window.
- The gap criterion is undecided. The best defect, 0.218(15), is within one error of the provisional
  scale 0.23; some leakage eigenvalues come out negative, so the leakage-sensitive figures need a
  bin-size check; and the omitted energy on the fine lattice is not measured here.

The masses (`blocked_masses.jl`, GEVP on levels 4 and 8, effective masses in fine units; coarse values
converted with a_c/a_f = 2.0015):

| fine t | α = 0.6, zero momentum | α = 0.3, zero momentum | coarse, starting at the same physical time (spans fine t → t + 2) |
|---|---|---|---|
| 0 | 0.768(13) | 0.794(13) | 0.540(7) (coarse t = 0) |
| 1 | 0.689(22) | 0.695(23) | |
| 2 | 0.656(40) | 0.666(42) | 0.486(17) (coarse t = 1) |
| 3 | 0.563(77) | 0.567(81) | |
| 4 | 0.529(133) | 0.534(139) | 0.457(44) (coarse t = 2) |

Published reference (Athenodorou and Teper, arXiv:2007.06422, Table 1): a m₀ = 0.987(9) at β = 5.6924
(10³ × 16), against our 0.973(35) at β = 5.7, and 0.7045(65) at β = 5.99, 0.6365(43) at β = 6.0625,
interpolating to about 0.645 at β = 6.054 in larger boxes than ours. In units of the string tension
their m₀/√σ rises from 2.47 at β = 5.6924 to 3.27 at β = 6.0625 (continuum 3.405(21)): the scalar dip of
the Wilson action. The rows of the table start at the same physical time but do not span the same
interval: a fine effective mass at t covers t → t + 1, the coarse one about t → t + 2 in fine units.
Over a matched interval, fine t = 2 → 4 against coarse t = 1 → 2, the central values are about 0.61 and
0.486 in fine units (α = 0.6, zero momentum); the error of that difference needs the jackknife over the
correlators and comes with the 3000-configuration analysis. The single fine effective masses fall from
0.656(40) at t = 2 to 0.56(8) and 0.53(13) at t = 3 and 4, consistent both with the benchmark and with
the coarse mass. A finite-spacing shift of the scalar mass, which
would account for most of the eigenvalue difference, is suggested but not established at 1000
configurations; resolving t = 3–5 is what the extension to 3000 is for. The decisive comparison is
over matched intervals: the fine mass over t = 2 → 4 (one mass from C(2)/C(4), about 0.61 at 1000
configurations) against the coarse mass over t = 1 → 2, and whether the difference persists one matched
interval later (`blocked_masses.jl` now prints this table). The printed overlaps are squared overlaps
with the leading variational state of the two-operator window, not with the true glueball, so
excited-state contamination remains a possible part of the early-time excess.

Scope. What is secure: no unitary alignment removes the observed one-step discrepancy within these
operator windows. The result concerns the tested blocked operators over one matched interval. A shift
of the scalar mass between the two spacings would obstruct a map that both pairs the two scalar ground
states and intertwines their evolution over matched times; pairing two rays is always possible, and
the constraint is that their decay rates agree. It leaves open how well other maps or other states
transfer. A lattice pair further into the scaling region would make a cleaner scaling study, while the
present pair can still test a finite-step inequality wherever its hypotheses hold. With the scale used
here, a factor 2 from β = 6.054 means β ≈ 6.55 (β = 6.4 is only 1.66, and a 32³ box there would be 21%
larger than the present one). The question the 3000-configuration run addresses: how much of the
one-step discrepancy persists in the long-time scalar decay?

## Full run at β = 6.054: 3000 configurations

The same chain extended to 3000 configurations (4.9 s per configuration on 12 threads, no projection
fallbacks). Bins of 10; the bin-size checks (20, 50) and the four-level window are recorded below when
available. The output file will be added to `lattice/results/`.

Masses (`blocked_masses.jl`, GEVP on levels 4 and 8, fine units). The fine effective mass has a plateau:
0.674(13), 0.658(24), 0.663(49) at t = 1, 2, 3 for α = 0.6 with zero momentum, and the same within errors
for α = 0.3 and for the fixed origin. It agrees with the published benchmark (about 0.645 at β = 6.054,
in larger boxes). Matched intervals, mass per fine step (fine t = 2j → 2j+2 against coarse j → j+1):

| j | fine, α = 0.6 zero | coarse | difference | α = 0.3 zero | α = 0.6 fixed | α = 0.3 fixed |
|---|---|---|---|---|---|---|
| 0 | 0.722(9) | 0.540(7) | 0.182(12) | 0.198(12) | 0.294(13) | 0.357(13) |
| 1 | 0.660(32) | 0.486(17) | 0.174(37) | 0.178(37) | 0.181(41) | 0.192(44) |
| 2 | 0.730(134) | 0.457(44) | 0.27(14) | 0.24(14) | 0.29(17) | 0.32(20) |

(The last three columns are differences.) The difference persists one matched interval later: for zero
momentum it is 0.182 at j = 0 and 0.174 at j = 1, and all four recipes agree at j = 1 (0.17–0.19, each
about 4.5 standard errors; they share configurations). The discrepancy therefore lies in the long-time
scalar decay, the mass itself, not only in short-time excited-state content. In the same physical units
the fine mass is 1.32 in coarse units against the coarse 0.973, about 36% heavier; with the published
string tension near β = 6.054 that is m₀/√σ ≈ 3.3, against about 2.4 on the coarse lattice. Both masses
agree with the published Wilson-action values, so the run is consistent with a finite-spacing shift of
the scalar mass between β = 5.7 and 6.054 (the scalar dip); attributing it to the bulk endpoint would
require varying the action.

Window {4, 8} (`window_transfer.jl --mref 0.973`):

| α, blocking | ‖K_f V − V K_c‖ | defect | total departure | eig(V*K_f V − K_c) | ε₊ | q floor | q needed (fine k) |
|---|---|---|---|---|---|---|---|
| 0, fixed origin | 0.301(6) | 0.333(9) | 0.183(7) | −0.050(6), −0.301(6) | 0 | 0.20(8) | 0.352(9) |
| 0, zero momentum | 0.252(6) | 0.286(10) | 0.165(8) | +0.010(6), −0.252(6) | 0.010(6) | 0.24(17) | 0.306(14) |
| 0.3, fixed origin | 0.208(6) | 0.256(10) | 0.156(9) | −0.003(6), −0.208(6) | 0 | 0.12(4) | 0.304(16) |
| 0.3, zero momentum | 0.173(6) | 0.209(12) | 0.129(11) | +0.040(6), −0.173(6) | 0.040(6) | 0.05(7) | 0.322(26) |
| 0.6, fixed origin | 0.193(6) | 0.243(10) | 0.150(9) | +0.008(6), −0.193(6) | 0.008(6) | 0.12(6) | 0.299(18) |
| 0.6, zero momentum | 0.166(6) | 0.203(13) | 0.126(11) | +0.045(6), −0.166(6) | 0.045(6) | 0.05(7) | 0.326(27) |

- Averaging over the blocking origins lowers the symmetric quantities at every α; at α = 0.6 the paired
  differences are 0.0273(25) for the in-window mismatch, 0.0396(53) for the defect and 0.0243(45) for the
  total departure.
- It also introduces a direction in which the fine evolution decays more slowly than the coarse one:
  ε₊ = 0.045(6) and 0.040(6) with zero momentum at α = 0.6 and 0.3, against 0.008(6) and 0 with the fixed
  origin. For a lower bound on the fine gap, the one-sided criterion, the fixed-origin map is the better
  one. This is a second axis on which the recipes rank in opposite order, after near-neutrality against
  near-invariance.
- Certificate. With the fine window measured directly, a fine gap at least the coarse reference
  (s = e^{−0.973}) needs q ≤ 0.30–0.33 for every recipe, above the positivity floors (0.05–0.24, noisy).
  With the fine window inferred from the coarse lattice no q suffices for any recipe: the coarse top
  0.339, plus ε₊, plus the leakage term leaves no room below s = 0.378.
- The best symmetric defect, 0.203(13), is 2.2 errors below the provisional scale 0.23 (within one error
  in the pilot); the symmetric defect is not the quantity a lower bound needs.

## One-sided bound

For a lower bound on the fine gap only one direction of disagreement matters: whether the fine
evolution could decay more slowly than the coarse one. Let 0 ≤ T ≤ 1 be the fine transfer matrix over
the matched interval with TΩ = Ω, P the window (orthogonal to Ω), Q the rest of Ω's complement, and

    k = λ_max(PTP),   ℓ² = ‖QTP‖² = λ_max(L_f),   q = ‖QTQ‖.

Writing a unit ψ ⊥ Ω as u + v (u = Pψ, v = Qψ) gives ⟨ψ, Tψ⟩ ≤ k‖u‖² + 2ℓ‖u‖‖v‖ + q‖v‖², hence

    λ₂(T) ≤ μ = (k + q)/2 + √(((k − q)/2)² + ℓ²),

the top eigenvalue of [k ℓ; ℓ q] (`one_sided_bound`). If q is well below k the leakage enters at
second order, μ ≈ k + ℓ²/(k − q); as q approaches k it enters at first order. To infer the fine window
from the coarse lattice only the one-sided part of the intertwining is needed: with
D = V*K_f V − K_c and ε₊ = max(0, λ_max(D)), k ≤ λ_max(K_c) + ε₊. Two constraints frame the omitted
top q:

- Positivity of T gives QTP = (QTQ)^{1/2} C (PTP)^{1/2} with ‖C‖ ≤ 1, so L_f ≤ q K_f and
  q ≥ λ_max(K_f^{-1/2} L_f K_f^{-1/2}) (`leakage_q_floor`). The data bound q from below; a gap
  certificate needs it bounded from above.
- μ ≤ s holds exactly when k < s and q ≤ s − ℓ²/(s − k) (`q_needed`).

Pilot (window {4, 8}, α = 0.6, zero momentum, central values): the eigenvalues of D are −0.1670 and
+0.0422, so most of the symmetric mismatch is in the harmless direction and ε₊ ≈ 0.042. The positivity
floor is q ≥ 0.0146/0.2331 ≈ 0.063 from the largest leakage, about 0.18 from the full matrix (noisy).
Certifying a fine gap equal to the coarse reference over the matched interval, s = e^{−0.973} ≈ 0.378,
needs q ≤ 0.277 with the directly measured k = 0.2331; with the coarse-derived k = 0.3392 + 0.0422 =
0.3814 > s no q suffices. The difference between the two shows what is lost when the fine window must be
inferred from the coarse lattice. In words, q ≤ 0.277 means that every normalized scalar state outside
the window retains at most that much autocorrelation over one interval. An omitted state could be
smooth and have the right symmetry; the question is whether the chosen operators see it. `window_transfer.jl` prints these quantities with jackknife errors
(`--mref 0.973` for this reference).

What an argument about the vacuum would have to supply. Not a smallness of the window's correlators:
a state orthogonal to every tested operator, with transfer eigenvalue arbitrarily close to 1, leaves
every window correlator at every time unchanged while q approaches 1, and a tiny overlap barely changes
a finite-time measurement. What is needed is a quantitative lower bound on how much of every
low-energy state, including linear combinations, the blocked operators capture; and since Q mixes
energy eigenstates, the operator to bound is Q e^{−τH} Q itself, not the exponential of a minimum mean
energy on Q. A result in the zero-momentum scalar sector controls that sector's gap; the full theory
needs the other sectors. Finite-range decompositions in rigorous renormalization exist even for the
massless Gaussian field, so they are a precedent for the form of such statements, not a substitute: the
missing step is the connection to this omitted Hilbert space. On the scalar dip: selective softening
of the scalar near the bulk endpoint of the fundamental–adjoint plane was found in mixed-action
simulations (Heller, hep-lat/9508009); the present two-ensemble run tests consistency of the mismatch
with a finite-spacing shift of the scalar mass, and attributing it to the endpoint would require
varying the action.

## The window map in the doubled space

Put the coarse and fine spaces side by side, H_c ⊕ H_f, with the doubled evolution S = S_c ⊕ S_f and,
for the window map W = U_f V U_c*, the polar reflection R = [0 W*; W 0] of the Universal Model (§3). On
the matched windows R is an involution; its +1 and −1 eigenspaces are the aligned graph Γ₊ = {(x, Wx)}
and the anti-aligned graph Γ₋ = {(x, −Wx)}. Multiplying the blocks gives

    ‖[S, R]‖ = ‖S_f W − W S_c‖,

so the defect of Proposition 6 (`window_defect`) measures how far Euclidean evolution fails to commute
with the operation that exchanges each matched coarse description with its fine partner. Two further
readings of the printed numbers (`graph_departure`; all checked on exact models in the tests):

- Normalization. For the raw correspondence A = V M induced by blocking, ‖M − I‖ is the largest change
  of length of a unit coarse window vector, and ‖M² − I‖ the largest absolute imbalance of squared
  lengths, |‖x‖² − ‖Ax‖²| for unit x: the largest absolute Krein square of the raw graph {(x, Ax)} with
  J = diag(I, −I), coarse first. Both vanish exactly when the raw map is isometric on the window; the
  polar graph is neutral by construction, and these numbers say how much normalization had to change.
- Evolution. One step of S carries the aligned graph partly into the anti-aligned graph, with amplitude
  ‖P₋SP₊‖ = ‖K_c − V*K_f V‖/2, half the in-window mismatch, and partly out of both windows, with
  ‖(I − P_c − P_f)SP₊‖ = λ_max((L_c + V*L_f V)/2)^{1/2}; the total departure ‖(I − P₊)SP₊‖ combines the
  two. A small in-window part with an appreciable outside part says the two descriptions agree on what
  the chosen operators retain while the operators miss part of the evolution.

Control: comparing the coarse ensemble with itself (levels 4, 8; V = I) gives 0 into the anti-aligned
graph and 0.122 both out of the windows and in total, the square root of the largest coarse leakage
eigenvalue 0.0149. So a total departure near 0.12 can coexist with perfect agreement; disagreement
inside the window is read from ‖K_f V − V K_c‖, whose pilot-size noise reference is 0.024 (0.012 for the
anti-aligned amplitude, its half). These are reference values, not calibrated thresholds for the fine
ensemble.

Two cautions. Centering an operator, O → O − ⟨O⟩, removes the vacuum component of the state it creates;
it is not the choice of the zero of energy, which the normalized correlators fix (they measure T/λ₀).
The relative phase of the graph is a free choice in the doubled description: nothing in these runs
prepares or detects it across lattices; what they test is whether the two descriptions stay paired
under evolution. At finite periodic time extent the ensemble correlators are thermal traces, and
reading them as vacuum matrix elements assumes vacuum dominance. Even one level gives an apparent
leakage from the backward term: C(t) ∝ r^t + r^{N_t − t} with am = 0.967 and N_t = 16 gives √L ≈ 0.001
at coarse lag 1, and the same at fine lag 2 with N_t = 32, the physical times being matched. This
illustrates the scale; it does not bound every thermal contribution. And pairing alone does not bound
the fine gap: with coarse energies 0, 1 and fine
energies 0, 1, ε, pairing the first two levels exactly gives a neutral graph that commutes with
evolution while the fine gap is ε. The omitted level is invisible to every comparison inside the
window, which is why Propositions 1 and 2 need the omitted energy Λ, the bottom of the fine Hamiltonian
on the complement of the window.

## Related work

Closest in structure: Faizal and Shabir, arXiv:2606.19362 (Fortschr. Phys. 74, 2026), organize a
construction of four-dimensional SU(N) Yang–Mills around block-spin maps: partial isometries V_k between
successive scales (a_{k+1} = b a_k, the coarse space embedded in the fine one, V_k Ω_{k+1} = Ω_k), an
interlacing identity T_{k+1} = V_k* T_k^b V_k − D_k + E_k, and summable defects. Their one-step estimate
(Eqs. 5.33–5.36, in physical units) bounds the gap of the coarser lattice from below by that of the
finer one. Corollary 3 runs the other way, a lower bound for ever finer lattices, and needs what theirs
does not: a lower bound on the energy of the fine states that blocking omits. Without it the inference
fails: with fine energies 0, ε, 1 and blocking retaining 0 and 1, the retained dynamics agree exactly
at matched times while the fine gap is ε. Two points about the statements. The hypotheses of Theorem 5.6
and Corollary D.5 give the telescoping bound Δ₀ − Σε_j; strict positivity needs the additional margin
Σε_j < Δ₀, which the proof of Corollary D.5 invokes, attributing it to a sufficiently large blocking
factor. Corollary 3 requires such a margin uniformly as the finest spacing shrinks. And the proof of
Theorem 11.1 uses the lower operator bound (11.3) as an upper bound in (11.5).

Parts of the measurement have precedents, two of them closest. Identifications between different
Hilbert spaces that control spectra are the subject of quasi-unitary equivalence (O. Post, Spectral
Analysis on Graph-like Spaces, Lecture Notes in Math. 2039, 2012; Post and Simmer, arXiv:2503.22611):
‖M − I‖ and the window defect play the roles of near-unitarity and approximate intertwining there, and
a detailed comparison is still to be made. The leakage matrix L = U*S²U − (U*SU)² of one lattice is,
with the symmetric normalization β₁ = γ₁ = C(0)^{1/2}, the residual Δ₂ = A₁(2) − α₁α₁ after the first
block-Lanczos iteration of Hackett and Wagman (Phys. Rev. D 112, 014514, 2025; arXiv:2412.04444,
Eq. 21), used there for spectroscopy on one lattice; here one such matrix enters on each side of a
blocking map fixed by the physics. Two-lattice Monte Carlo renormalization group matching, a method of the 1980s, has been
applied to SU(3) gauge theories (A. Hasenfratz, arXiv:0907.0919). Matching two lattice spacings through
glueball correlations was itself considered then: glueball correlation functions or the heavy-quark
potential could in principle be used, but their long-wavelength properties are statistically hard, so
blocked loop expectation values were matched instead; and no blocking kernel is expected to lead to
exactly the Wilson action at finite spacing, since that would make the Wilson action free of
finite-spacing errors (P. B. Mackenzie, Monte Carlo renormalization group calculations for SU(2)
lattice gauge theory, FERMILAB-Conf-84/48-T, 1984). The question was seen; the
aim then was to tune the coarse action, with observables chosen for that. The intertwining relation
S_f W = W S_c has the form of the lumpability condition QC = Cq for Markov chains, and Ward et al.
(Commun. Phys., 2022, doi:10.1038/s42005-022-00834-1) measure the failure of approximate lumping by
‖QC − Cq‖, the algebraic form of the defect; there the coarse generator minimizes its Frobenius norm
(their Theorem 2.1), whereas here the operator norm is measured, with the coarse dynamics taken from an
independently simulated lattice. Lang, Liegener and Thiemann (Hamiltonian renormalisation I–IX, arXiv:1711.05685 and
sequels) derive the compression of the semigroup to the coarse space that blocking induces through
Osterwalder–Schrader reconstruction, and when the embedded space is invariant; their work is analytic,
for free and solvable models. Isometric scaling maps between lattice Hilbert spaces with an inductive limit: Morinelli,
Morsella, Stottmeister and Tanimoto (Commun. Math. Phys. 2021), Stottmeister et al. (Phys. Rev. Lett.
127, 230601, 2021), Brothier and Stottmeister for 1+1-dimensional Yang–Mills (arXiv:1907.05549).
Effective transfer matrices on block variables estimated by Monte Carlo, with an overlap measure of
quality, on one lattice: Hasenbusch, Pinn and Rummukainen (Phys. Rev. D 49, 495, 1994). Standard
practice also compares blocked observables through fixed-point actions, blocked and smeared glueball
operators, and gradient flow as a continuous blocking (e.g. Carosso, Hasenfratz and Neil,
arXiv:1806.01385). Rigorous background: Bałaban's ultraviolet stability for lattice Yang–Mills
(expounded by Dimock) and the strong-coupling gap of Osterwalder and Seiler.

What we put forward is a test of whether the correspondence between two lattice descriptions preserves
what happens to the states it pairs, a physics question that stands whether or not a proof of the gap
follows. Evolving a coarse state one step and then translating it, or translating it first and evolving
the fine state over the same physical interval, the defect is how far apart the two end up. The
measurement is the combination: blocking fixes which coarse operator corresponds to which fine
operator; the polar map aligns their normalizations; correlators measured independently at the two
spacings then determine the full intertwining defect of that window map, including evolution out of
either window, and show how much of the disagreement is already visible inside the paired windows. We
found precedents for each ingredient, but not this measurement for SU(3), in the sources checked (a web
search and the closest papers, not a systematic citation search). The physics question it answers: when
two lattice spacings give compatible glueball masses, do the corresponding operators also evolve alike?
Two operators can share the leading decay while weighting higher excitations differently, and the
measurement exposes that difference; Corollary 3 says why its size matters. Whether a disagreement can
weaken a lower bound on the gap depends on its direction (the one-sided bound). The clearest payoff of
the polar formalism so far is the pilot's opposite ranking at α = 0.6: the fixed grid preserves the
lengths of the paired states better, the origin-averaged grid makes their evolution agree better.
Judged by how faithfully it prepares the state at the initial instant, the map would have been chosen
wrongly for tracking the state one step later; separating near-neutrality of the graph from its
near-invariance is what told us to measure the two apart.

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

Newest: the full run ("Full run at β = 6.054"): the persistence of the mass difference over
matched intervals, and the one-sided reversal between the recipes. Before that: the one-sided bound ("One-sided bound"): the 2 × 2 lemma, the positivity floor on q, the q
needed for a reference gap, and what a vacuum argument would have to supply. Before that: the pilot results at β = 6.054 ("Pilot at β = 6.054"): the spectral reading of the mismatch,
the mass comparison and its interpretation, and the scope statement. Before that: the doubled-space
reading and `graph_departure` (the identity ‖[S, R]‖ = ‖S_fW − WS_c‖ and
the three departure formulas), and the related-work section, in particular the reading of
arXiv:2606.19362.

Reliability patch (previous round): the checkpoint record and resume checks, `read_measurements`, the
branch-cut band of `su3_project` and its tests, `lattice/test/runner_checks.jl`, and the paired
fixed − zero differences in `window_transfer.jl`, all described under "Protecting the measurements".

Previous round:

1. Proposition 6 and Remark 7: the bound and the equality for the defect itself, the alignment induced
   by blocking as the polar part of the operator-state correspondence, centered excitations, and the
   fixed-origin and zero-momentum variants of the blocked operators.
2. `block_spatial` (now with an origin) and `su3_project` (with a covariant fallback, the straight path,
   for numerically singular path sums), `window_defect`, and the dimension and row-width checks of
   `window_transfer.jl`.
3. The pilot and budget above: whether the fine run, as designed, can answer the question.
4. `examples/blocked_glueball.jl` and `examples/window_transfer.jl`: the two-ensemble jackknife and the
   use of fine lag 2 against coarse lag 1.

Earlier:

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
