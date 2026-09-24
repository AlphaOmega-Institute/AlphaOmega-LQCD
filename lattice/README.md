# LatticeYM

A small, tested engine for pure SU(3) Yang–Mills theory on the lattice, in Julia with no dependencies outside the standard library.

- Wilson plaquette action on a periodic hypercubic lattice in any number of dimensions
- Optional 't Hooft twisted boundary conditions
- Cabibbo–Marinari heatbath (Kennedy–Pendleton) with overrelaxation, and multi-hit Metropolis
- Plaquette (total, per plane), Polyakov loop, R×T Wilson loops
- Error bars that account for autocorrelation (Madras–Sokal windowing)

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

The statistical tests allow 4.5 standard errors, with errors that include autocorrelation.

To check the tests would catch real mistakes, each of these bugs was planted by hand, and each made the suite fail:

| Planted bug | Failing checks |
|---|---|
| Twist phase the same for both orientations (the rule in the older `src/julia/viz/check_pi_4.jl`) | 18: staple vs action, corner vs bulk, overrelaxation |
| Twist left out of the update, kept in the measurement | 31: same tests |
| Heatbath using half the correct coupling | 7: every exact-plaquette comparison |
| Metropolis proposals that are not symmetric | 4: the exact 2D values and the heatbath cross-check |

## Results

Reproduce with `lattice/examples/beta_scan.jl`; a table of reference runs will be added here.

## Not yet included

- Gauge groups other than SU(3), and the Twisted Eguchi–Kawai (single-site, large-N) model
- Gradient flow, smearing, topological charge, scale setting
- Performance work beyond avoiding allocations: about 1 µs per link update on one core (heatbath plus staple)
