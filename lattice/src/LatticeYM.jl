"""
    LatticeYM

Pure SU(3) Yang–Mills theory on a periodic hypercubic lattice in any number of
dimensions, with the Wilson plaquette action and optional 't Hooft twisted
boundary conditions.

Every piece is checked against exact results in `test/runtests.jl`; see the
package README for what is and is not validated.
"""
module LatticeYM

using LinearAlgebra
using Random

export Mat3, retr, reunitarize, haar_su3, su3_deviation
export Lattice, GaugeField, cold_start, hot_start, nsites
export staple, plaquette, plaquette_retr, action, gauge_transform!
export mean_plaquette, plane_plaquettes, corner_and_bulk, polyakov_loop, wilson_loop
export Metropolis, Heatbath, Overrelaxation, sweep!
export mean_err, single_plaquette_average
export stout_spatial, slice_plaquettes, glueball_operators, wilson_loops_rt
export bin_means, jackknife, correlator_matrix, correlator_products, correlator_from_means
export cosh_mass, effective_masses, gevp, transfer_residual
export su3_project, block_spatial
export window_compression, induced_alignment, window_defect_bound, window_defect, ground_state_fractions

include("su3.jl")
include("geometry.jl")
include("gaugefield.jl")
include("updates.jl")
include("observables.jl")
include("analysis.jl")
include("exact.jl")
include("smearing.jl")
include("spectroscopy.jl")

end
