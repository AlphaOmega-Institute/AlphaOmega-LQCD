# Reference results

Curated datasets behind results quoted in the documentation, kept so that the analysis can be
reproduced without rerunning the simulation. Routine simulation output belongs in `data/`, which git
ignores.

| File | Contents |
|---|---|
| `glueball_b570_L8T16.txt.gz` | Output of `examples/glueball.jl` at β = 5.7 on 8³ × 16: 3000 configurations, each with the mean plaquette, the 0⁺⁺ operators at stout levels 0, 4, 8, 16, 24, 32 (ρ = 0.1) on all 16 time slices, and Wilson loops R ≤ 4, T ≤ 6 at level 16. The command and seed are in the header line. |
| `glueball_b570_L8T16_analysis.txt` | `examples/glueball_analysis.jl` on that file with `--bin 10 --ops 2,3,4,5`. |
| `fine_b6054.txt.gz` | Output of `examples/blocked_glueball.jl` at β = 6.054 on 16³ × 32, blocked once to 8³ with α = 0, 0.3, 0.6 at all eight origins: 3000 configurations (4 sweeps apart after 500 thermalization sweeps, seed 6054), each with the fine mean plaquette and, for each α, the fixed-origin and zero-momentum 0⁺⁺ operators at stout levels 4, 8, 16, 24 (ρ = 0.1) on all 32 time slices. One chain, run at commit 4ec40b3 on 12 threads of an Apple M4 Max, checkpointed and resumed at configuration 1000; no projection fallbacks. The header line records the options. |
| `fine_b6054_masses.txt` | `examples/blocked_masses.jl` on that file with `--coarse` the β = 5.7 file: effective masses, matched intervals and correlated fits. |
| `fine_b6054_window48_bin10.txt`, `_bin20`, `_bin50` | `examples/window_transfer.jl` with the β = 5.7 file as coarse, `--window 4,8 --mref 0.973`, bins of 10, 20 and 50. |
| `fine_b6054_window4to24_bin10.txt` | The same with `--window 4,8,16,24`, bins of 10. |

The analysis outputs were first produced on the machine that ran the simulation (commit f3ade0c);
rerun from the uploaded data on a different machine they were reproduced byte for byte. The files here
come from the current scripts, which add the correlated fits and the paired difference of ε₊.
