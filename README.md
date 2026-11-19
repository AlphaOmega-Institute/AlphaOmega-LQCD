# AlphaOmega-LQCD

[![Julia](https://img.shields.io/badge/Julia-1.9+-9558B2?logo=julia&logoColor=white)](https://julialang.org)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17844109.svg)](https://doi.org/10.5281/zenodo.17844109)

Lattice gauge theory code for the A|Ω⟩ research program.

The repository currently contains one package, [**LatticeYM**](lattice/README.md): pure SU(3) Yang–Mills theory on the lattice, written in Julia with no dependencies outside the standard library.

- Wilson plaquette action on a periodic lattice in any number of dimensions, with optional 't Hooft twisted boundary conditions
- Heatbath with overrelaxation, and Metropolis
- Plaquettes, Polyakov loops, Wilson loops, and error bars that account for autocorrelation
- A test suite that checks the code against exact results: the 2D plaquette, 4D strong coupling, gauge invariance, a full recomputation of the action, and the published β = 6.0 plaquette

This is standard lattice gauge theory, meant as a tested foundation for later work.

## Quick start

```sh
julia --project=lattice -t auto lattice/test/runtests.jl                      # about 1 minute
julia --project=lattice -t auto lattice/examples/beta_scan.jl --L 8 5.7 6.0   # plaquette at given β
```

See [lattice/README.md](lattice/README.md) for conventions, what each test checks, and reference results.
