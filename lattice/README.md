# LatticeYM

A small, tested engine for pure SU(3) Yang–Mills theory on the lattice, in Julia with no dependencies outside the standard library.

- Wilson plaquette action on a periodic hypercubic lattice in any number of dimensions
- Optional 't Hooft twisted boundary conditions
- Cabibbo–Marinari heatbath (Kennedy–Pendleton) with overrelaxation, and multi-hit Metropolis
- Plaquette (total, per plane), Polyakov loop, R×T Wilson loops
- Error bars that account for autocorrelation (Madras–Sokal windowing)
- Glueball spectroscopy: spatial stout smearing, 0⁺⁺ operators per time slice at several smearing levels, smeared R×T Wilson loops for the static potential, and a correlator analysis (binning, jackknife, periodic effective masses, generalized eigenvalue problem)
- Spatial blocking with a gauge-covariant SU(3) projection, and operator-window quantities for comparing a coarse lattice with a blocked fine one

This is standard lattice gauge theory, not new physics. It is meant as a foundation whose results can be trusted, because every piece is checked against something exact.

## Running

```sh
julia --project=lattice -t auto lattice/test/runtests.jl                       # test suite, about 1 minute
julia --project=lattice -t auto lattice/examples/beta_scan.jl --L 8 5.7 6.0    # plaquette at given β
```

```julia
using LatticeYM, Random
rng = Xoshiro(1)
lat = Lattice((8, 8, 8, 8); twist = ((1, 2) => 1, (3, 4) => 1))   # twist is optional
g = hot_start(lat, rng)
for _ in 1:500
    sweep!(g, 6.0, Heatbath(), rng)
    sweep!(g, 6.0, Overrelaxation(), rng)
end
mean_plaquette(g)
```

### Glueball spectroscopy

```sh
julia --project=lattice -t auto lattice/examples/glueball.jl --beta 5.7 --L 8 --T 16 --confs 3000 \
      --levels 0,4,8,16,24,32 --wlevel 16 > data.txt
julia --project=lattice lattice/examples/glueball_analysis.jl data.txt --bin 10 --ops 2,3,4,5
```

To compare two lattice spacings, block a fine lattice and measure the same operators, then compare
with the coarse data (`docs/notes/spectral_transfer.tex`, Proposition 6):

```sh
julia --project=lattice -t auto lattice/examples/blocked_glueball.jl --checkpoint fine.ckpt > fine.txt
julia --project=lattice lattice/examples/window_transfer.jl data.txt fine.txt --window 4,8
```

The first script writes, per configuration, the mean plaquette, the 0⁺⁺ slice operators at each smearing
level and smeared R×T Wilson loops. The second reports effective masses per level, the variational
(GEVP) energies, ground-state fractions by two estimators, the defect bound of
`docs/notes/spectral_transfer.tex`, the transfer-matrix residual of each trial state, the static potential
and the string tension.

## Conventions

- `U[μ, x]` is the link from site x to x + μ̂.
- P_μν(x) = U_μ(x) U_ν(x+μ̂) U_μ(x+ν̂)† U_ν(x)†.
- S = β Σ_{x, μ<ν} (1 − Re tr(z P_μν(x)) / 3), so the action density is ⟨S⟩ = 1 − ⟨P⟩.
- Twist `(μ, ν) => n`: in every (μ, ν) plane-slice, the single plaquette at x_μ = L_μ, x_ν = L_ν carries z = exp(2πi n/3), and the opposite orientation carries conj(z). All other plaquettes have z = 1.

## What the tests check

| Test | Compared against |
|---|---|
| SU(3) algebra, reunitarization, SU(2) embeddings | Dense-matrix arithmetic; unitarity and det = 1 to 1e-13 |
| Haar sampler | ⟨tr U⟩ = 0 and ⟨\|tr U\|²⟩ = 1 |
| Heatbath a₀ sampler | Exact mean of √(1−a²) e^{αa}, for α from 0 to 12 |
| Staple | Change in the full action after replacing one link, with and without twist, including links at the twisted corner |
| Gauge invariance | Action, plaquettes, Polyakov and Wilson loops unchanged under a random gauge transformation, with and without twist |
| Overrelaxation | Leaves the total action unchanged |
| 2D plaquette at β = 1, 4, 8 | Exact result from the Weyl integral, for both heatbath and Metropolis |
| 2D with twist | Same exact value, both at the twisted plaquette and away from it |
| 4D at β = 1 | Strong-coupling result (corrections of order 10⁻⁶) |
| 4D at β = 5.7 | Heatbath and Metropolis agree |
| 4D with twist | Twisted corner plaquettes average the same as the rest of their plane, as they must, because a change of variables moves the twist |
| Reproducibility | Same seed gives identical results |
| Matrix exponential | `expm3` against the dense exponential to 1e-12, for norms from 0.01 to 12; result in SU(3) |
| Stout smearing | Links stay in SU(3); time links unchanged; unit and pure-gauge fields stay so; ρ = 0 is the identity; smearing commutes with random gauge transformations to 1e-12; spatial plaquettes increase |
| Glueball operators | Slice average equals the spatial plane plaquettes; each smearing level matches direct smearing; gauge invariant |
| Smeared Wilson loops | `wilson_loops_rt` equals `wilson_loop` averaged over spatial directions; W(1,1) equals the space-time plaquettes; gauge invariant |
| SU(3) projection | `su3_project` lands in SU(3), fixes SU(3) matrices and their positive multiples, and satisfies su3_project(G A H†) = G su3_project(A) H† to 1e-11; row Gram–Schmidt fails that test |
| Spatial blocking | Halved spatial extents; links in SU(3); α = 0 gives the product of two links; the unit field stays unit; blocking commutes with gauge transformations restricted to coarse sites; operators on blocked fields are gauge invariant |
| Blocking origins | Each origin maps coarse sites to the right fine sites and is gauge covariant; a one-site translation changes the fixed-origin operator but leaves the eight-origin average unchanged; singular path sums fall back to the straight path |
| Operator windows | `window_defect` equals ‖S_f W − W S_c‖ on exact models; `window_compression` reproduces U*SU and U*S(1 − P)SU on exact models, for lag 1 and 2; `induced_alignment` gives the polar factorization; the Proposition 6 bound holds; `ground_state_fractions` is exact for two-state correlators at t₀ = 0, 1, 2, with and without periodic images |
| Correlator analysis | Effective masses of exact cosh correlators; GEVP recovers both energies of an exact two-state correlator at several (t₀, t); `correlator_matrix` against a direct implementation; binned products reproduce it; the jackknife error of a mean equals the standard error; `transfer_residual` equals ‖(I − P_ψ)Tψ‖ for random positive transfer matrices |

The statistical tests allow 4.5 standard errors, with errors that include autocorrelation.

To check the tests would catch real mistakes, each of these bugs was planted by hand, and each made the suite fail:

| Planted bug | Failing checks |
|---|---|
| Twist phase the same for both orientations | 18: staple vs action, corner vs bulk, overrelaxation |
| Twist left out of the update, kept in the measurement | 31: same tests |
| Heatbath using half the correct coupling | 7: every exact-plaquette comparison |
| Metropolis proposals that are not symmetric | 4: the exact 2D values and the heatbath cross-check |

## Results

Produced with `examples/beta_scan.jl`: hot start, seed 1, and 200 thermalization plus 2000 measured sweeps unless noted. One sweep is one heatbath pass and three overrelaxation passes.

| β | Lattice | ⟨S⟩ = 1 − ⟨P⟩ |
|---|---|---|
| 3.10 | 8⁴ | 0.786983(41) |
| 3.12 | 8⁴ | 0.785252(41) |
| 3.14 | 8⁴ | 0.783620(43) |
| 3.90 | 8⁴ | 0.718664(53) |
| 3.95 | 8⁴ | 0.714142(49) |
| 3.98 | 8⁴ | 0.711328(51) |
| 3.98 | 8⁴, twisted | 0.711431(55) |
| 4.00 | 8⁴ | 0.709587(54) |

The twist leaves the bulk action unchanged: twisted and periodic runs at β = 3.98 differ by 0.00010(8).

Check against the literature at β = 6.0, where the commonly quoted large-volume plaquette is ≈ 0.59368:

| Lattice | Sweeps (therm + measured) | ⟨P⟩ |
|---|---|---|
| 8⁴ | 200 + 2000 | 0.594374(88) |
| 12⁴ | 200 + 1000 | 0.593673(56) |
| 16⁴ | 150 + 400 | 0.593653(55) |

The 8⁴ box is too small at this coupling and sits about 0.0007 high. From 12⁴ upward the result agrees with the literature.

### Glueball spectroscopy at β = 5.7

8³ × 16, 3000 configurations, data in `results/glueball_b570_L8T16.txt.gz` (with the analysis output next
to it). String tension from smeared Wilson loops: σa² = 0.1677(7) from a Cornell fit with the Coulomb
coefficient fixed at π/12, and 0.1576(24) with it free; the spread between the two fits is larger than
either error. Lucini and Teper (hep-lat/0103027) find 0.1569(18) at this β from flux loops. 0⁺⁺ mass
a·m = 0.967(33) from a variational analysis, against their 0.999(25), so m/√σ = 2.36(8). This is
below the continuum value of about 3.5, as expected at this coupling: the scalar glueball has a known dip
near β = 5.7. Ground-state fractions per smearing level and the transfer-matrix residuals are tabulated
in `../docs/notes/README.md`.

## Not yet included

- Gauge groups other than SU(3), and the Twisted Eguchi–Kawai (single-site, large-N) model
- Gradient flow, topological charge, scale setting beyond the string tension
- Glueball channels other than A₁⁺⁺, and operators other than smeared plaquettes
- Performance work beyond avoiding allocations: about 1 µs per link update on one core (heatbath plus staple)
