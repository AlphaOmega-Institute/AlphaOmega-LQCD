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
| `fine_b6054_window48_mref1311_bin10.txt` | The first pair with `--window 4,8 --mref 1.311`, the reference mass of the second pair in β = 5.7 units, for the comparison at a common interval. |
| `coarse_b6054.txt.gz` | Output of `examples/glueball.jl` at β = 6.054 on 16³ × 32, unblocked: 3000 configurations (4 sweeps apart after 500 thermalization sweeps, seed 60540), the 0⁺⁺ operators at stout levels 4, 8, 16, 24 on all 32 time slices and Wilson loops R ≤ 6, T ≤ 8 at level 16. The second pair's coarse ensemble. |
| `fine_b6546.txt.gz` | Output of `examples/blocked_glueball.jl` at β = 6.546 on 32³ × 64, blocked once to 16³ with α = 0.3, 0.6 at all eight origins: 1000 configurations (8 sweeps apart after 1000 thermalization sweeps, seed 6546), operators at stout levels 4, 8, 16, 24 on all 64 time slices. Run at commit 675f3ee on 12 threads of an Apple M4 Max; no projection fallbacks. |
| `fine_b6546_window48_bin10.txt`, `_bin20`, `_bin50` | `examples/window_transfer.jl` with `coarse_b6054` as coarse, `--window 4,8 --mref 0.655`. |
| `fine_b6546_window48_lag2_bin10.txt` | The same with `--lag 2` (the first pair's physical interval). |
| `fine_b6546_window4to24_bin10.txt` | The same with `--window 4,8,16,24`, bins of 10. |
| `fine_b6546_masses.txt`, `_bin20` | `examples/blocked_masses.jl` with `--coarse coarse_b6054 --ratio 1.9997 --franges 4:16,8:16,4:2:16,2:8 --cranges 2:8,4:8,1:4`. |
| `fine_b6546_masses_bin50.txt` | The same at bins of 50 without the short ranges (`--franges 4:16,8:16,4:2:16 --cranges 2:8,4:8`). |
| `noise_b6054_3000v1000_window48.txt`, `_window4to24` | `examples/window_noise_floor.jl` on `coarse_b6054` with `--size 3000 --size2 1000`: the noise reference for the second pair. |
| `noise_b6054_3000v1000_window48_lag2.txt` | The same at `--lag 2`, the noise reference at the first pair's physical interval. |

The analysis outputs were first produced on the machine that ran the simulation (commit f3ade0c);
rerun from the uploaded data on a different machine they were reproduced byte for byte. The files here
come from the current scripts, which add the correlated fits, the paired difference of ε₊ and, after the
second pair's data were seen, how far each maximization lies above its value at D = 0.
