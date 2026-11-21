# AlphaOmega-LQCD

[![Julia](https://img.shields.io/badge/Julia-1.9+-9558B2?logo=julia&logoColor=white)](https://julialang.org)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17844109.svg)](https://doi.org/10.5281/zenodo.17844109)

Lattice gauge theory code for the A|Ω⟩ research program.

The repository currently contains one package, [**LatticeYM**](lattice/README.md): pure SU(3) Yang–Mills theory on the lattice, written in Julia with no dependencies outside the standard library.

- Wilson plaquette action on a periodic lattice in any number of dimensions, with optional 't Hooft twisted boundary conditions
- Heatbath with overrelaxation, and Metropolis
- Plaquettes, Polyakov loops, Wilson loops, and error bars that account for autocorrelation
- Glueball spectroscopy: stout smearing, 0⁺⁺ correlators, the static potential and the string tension
- A test suite that checks the code against exact results: the 2D plaquette, 4D strong coupling, gauge invariance, a full recomputation of the action, and the published β = 6.0 plaquette

This is standard lattice gauge theory, meant as a tested foundation for later work.

Also in this repository:

- [`hamiltonian/`](hamiltonian/README.md): exact SU(2) Hamiltonian models on one and two plaquettes (Python)
- [`docs/notes/`](docs/notes/README.md): working notes, currently on spectral transfer along block-spin maps

## Quick start

```sh
julia --project=lattice -t auto lattice/test/runtests.jl                      # about 1 minute
julia --project=lattice -t auto lattice/examples/beta_scan.jl --L 8 5.7 6.0   # plaquette at given β
```

See [lattice/README.md](lattice/README.md) for conventions, what each test checks, and reference results.

## License

Copyright © 2026 Alpha Omega Lattice Corporation, doing business as AlphaOmega Institute. The Julia code is licensed under the PolyForm Noncommercial License 1.0.0, and the documentation under CC BY-NC 4.0. The logo is not licensed. See [LICENSE.md](LICENSE.md) for details. Versions whose README showed an MIT badge, and the v0.1.0 archive on Zenodo, keep their earlier terms. For commercial uses that require permission beyond the applicable license, contact Alpha Omega Lattice Corporation.
