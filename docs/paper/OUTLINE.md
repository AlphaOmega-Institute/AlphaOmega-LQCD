# Paper outline (working)

Status: outline only. Numbers refer to the review guide, `docs/notes/README.md`, and to the files in
`lattice/results/`. The second lattice pair (β = 6.054 against 6.546) is planned; its prediction is
written down in the guide before the run.

## Working titles

- Does a block-spin map preserve the dynamics? A two-spacing test with SU(3) glueball operators
- Testing coarse-to-fine maps between lattice spacings through operator windows

## The claim, in one paragraph

Blocking fixes which coarse operator a fine operator represents. Evolving a coarse state over one step
and then translating it, or translating it first and evolving the fine state over the same physical
time, should agree if the map preserves the dynamics. We measure the largest disagreement over a window
of operators, exactly, from correlators of two independently simulated lattices. Three quantities have
distinct meanings, though their norms do not add: the minimum mismatch over all alignments (a floor set
by the spectra), the mismatch of the alignment blocking actually selects, and evolution out of the
window, which enters the full defect through an exact block-matrix formula. On
β = 5.7 against 6.054 the physical pairing reveals what the spectra hide: both fine window eigenvalues
lie below the coarse ones, yet the blocking-fixed pairing rotates the leading direction by about 25°
and creates a direction in which the fine evolution decays more slowly. We state the one-sided criterion
that a lower bound on a gap needs, and what it leaves for an independent argument: an upper bound on
the slow states the chosen operators do not see.

## 1. Introduction

- The two orders of "block and evolve", and the measured consequence up front: both sorted fine
  eigenvalues lie below their coarse partners, yet the physical pairing gives a direction with slower
  fine decay (a 25° rotation). Why comparing masses leaves something out.
- The question: does a map between two lattice descriptions preserve what happens to the states it
  pairs? A physics question in its own right, and the hypothesis of renormalization-group arguments for
  a mass gap (gap transfer along a chain of maps).
- What is usually compared between spacings: masses and other numbers; Monte Carlo renormalization
  group matching tunes the coarse action (Mackenzie 1984; Hasenfratz). What is new here: a map fixed in
  advance, and a comparison of dynamics, not of numbers.
- Results in brief; scope (one channel, finite volume, Wilson action; no proof of anything about the
  continuum).

## 2. Framework

- Window maps: partial isometries between window spaces spanned by operator states; compression K and
  leakage L from correlators at lags 0, ℓ, 2ℓ (the one-lattice L is the first block-Lanczos residual of
  Hackett and Wagman).
- Alignment fixed by blocking: polar part V of G_f^{1/2} G_c^{-1/2}; normalization M.
- Exact defect of the window map (Proposition 6), and the bound; spectral floor by Weyl's inequality.
- Doubled-space reading, briefly (derivation in an appendix): polar reflection R,
  ‖[S, R]‖ = ‖S_f W − W S_c‖; near-neutrality against near-invariance.
- One-sided criterion: the 2 × 2 lemma, positivity floor L ≤ qK, the q needed for a reference gap,
  the matrix form k ≤ λ_max(K_c + D₊).
- Gap transfer along a chain: direction, omitted energy, uniform margin; the operator form carried from
  rung to rung (proof in an appendix); relation to the related work (Faizal and Shabir; quasi-unitary
  equivalence, Post and Simmer).

## 3. Lattice setup

- Wilson action, heatbath and overrelaxation; stout smearing of spatial links; 0⁺⁺ slice operators.
- Gauge-covariant spatial blocking: path sums with polar SU(3) projection, fallbacks for singular sums
  and at the determinant branch cut; eight origins and the zero-momentum average.
- Ensembles (table): β, volume, configurations, separation, seeds; scale from the Necco–Sommer
  interpolation; matched physical volume.
- Analysis: binning, two-ensemble jackknife, correlated cosh fits with jackknife covariance; checks
  (bin size, window size, byte-for-byte reproduction on a second machine).

## 4. Results: β = 5.7 against β = 6.054

- Masses: effective masses, matched intervals, correlated fits; the fine-minus-coarse difference
  0.171(29) persists past the first interval; agreement with published values; the scalar dip accounts
  for the leading eigenvalue separation (0.110 predicted, 0.103 observed).
- Window quantities by recipe (table): mismatch, spectral floor, defect, departure, ε₊.
- What the pairing shows beyond the spectra: the 25° rotation and ε₊ = 0.045(6).
- Two reversals between recipes: length preservation against evolution tracking; symmetric agreement
  against the one-sided error; with leakage, the certificate ranking is open.
- One-sided certificate: q needed with the fine window measured directly (0.30–0.35), with the matrix
  form from the coarse lattice (0.146(125) for the primary recipe), positivity floors.

## 5. Results: β = 6.054 against β = 6.546

- The prediction and analysis plan, as written before the fine data; the results at one coarse step
  and at the first pair's physical interval; whether the angle, ε₊ and the leakage shrink with the
  spectral separation.

## 6. Discussion

- What the method measures that masses do not.
- For renormalization-group arguments: the direction of the inequality, operator inequalities rather
  than norms, the omitted energy as the hypothesis simulation cannot supply, and what an independent
  argument would have to show (a quantitative lower bound on how much of every low-energy state the
  blocked operators capture).
- Limitations: one channel, finite volume, the Wilson action, finite-time windows, statistics.

## 7. Conclusion

## Appendices

- Proofs: Proposition 6 (equality), the one-sided lemma, the positivity floor, the matrix-form bound,
  the operator-valued step for a chain.
- The doubled-space derivation.
- Code and data availability: repository, data files, commands; reproduction.
- Disclosure of the use of AI assistance, following the venue's policy.

## Figures and tables

1. Schematic of the two orders of evolving and translating.
2. Effective masses of both lattices at matched physical times, with published values.
3. Table of ensembles.
4. Table of window quantities by recipe.
5. Leading directions under spectral and blocking-fixed pairing (the rotation).
6. μ(q) for the direct and coarse-derived window tops, with the positivity floor and the target.
7. Second pair: prediction against result.

## Before drafting prose

- The second lattice pair.
- A systematic citation search (INSPIRE: works citing Hasenbusch–Pinn–Rummukainen, Mackenzie 1984,
  Hackett–Wagman, the Hamiltonian renormalisation series).
- Venue (a lattice-methods paper in a physics journal, or conference proceedings first), arXiv
  endorsement, authorship and affiliation, the venue's policy on AI assistance.
