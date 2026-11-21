# Reference results

Curated datasets behind results quoted in the documentation, kept so that the analysis can be
reproduced without rerunning the simulation. Routine simulation output belongs in `data/`, which git
ignores.

| File | Contents |
|---|---|
| `glueball_b570_L8T16.txt.gz` | Output of `examples/glueball.jl` at β = 5.7 on 8³ × 16: 3000 configurations, each with the mean plaquette, the 0⁺⁺ operators at stout levels 0, 4, 8, 16, 24, 32 (ρ = 0.1) on all 16 time slices, and Wilson loops R ≤ 4, T ≤ 6 at level 16. The command and seed are in the header line. |
| `glueball_b570_L8T16_analysis.txt` | `examples/glueball_analysis.jl` on that file with `--bin 10 --ops 2,3,4,5`. |
