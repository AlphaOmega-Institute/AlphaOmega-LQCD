# Measure 0⁺⁺ glueball operators at several smearing levels, and smeared R×T
# Wilson loops, on configurations of SU(3) Yang–Mills theory.
#
#   julia --project=lattice -t auto lattice/examples/glueball.jl [options] > data.txt
#
# Options: --beta 5.7        coupling
#          --L 8 --T 16      spatial and temporal extents (L³ × T)
#          --therm 300       thermalization sweeps
#          --confs 1000      measured configurations
#          --skip 4          sweeps between measurements
#          --rho 0.1         stout parameter (spatial links only)
#          --levels 0,5,10,20,30,45   smearing levels of the glueball operators
#          --wlevel 20       smearing level of the Wilson loops
#          --rmax 4 --tmax 6 largest Wilson loop
#          --seed 1
# One sweep is one heatbath pass and three overrelaxation passes. Output: a
# header, then one line per configuration with the mean plaquette, the operators
# O[level, t] (level-major), and W[R, T] (R-major).

using LatticeYM, Random, Printf

function main(args)
    opt = Dict("beta" => "5.7", "L" => "8", "T" => "16", "therm" => "300", "confs" => "1000",
               "skip" => "4", "rho" => "0.1", "levels" => "0,5,10,20,30,45", "wlevel" => "20",
               "rmax" => "4", "tmax" => "6", "seed" => "1")
    i = 1
    while i <= length(args)
        key = lstrip(args[i], '-')
        haskey(opt, key) || error("unknown option $(args[i])")
        opt[key] = args[i+1]
        i += 2
    end
    β, ρ = parse(Float64, opt["beta"]), parse(Float64, opt["rho"])
    L, T = parse(Int, opt["L"]), parse(Int, opt["T"])
    ntherm, nconf, nskip = parse.(Int, (opt["therm"], opt["confs"], opt["skip"]))
    levels = parse.(Int, split(opt["levels"], ","))
    wlevel, rmax, tmax, seed = parse.(Int, (opt["wlevel"], opt["rmax"], opt["tmax"], opt["seed"]))
    all_levels = sort(unique(vcat(levels, wlevel)))

    lat = Lattice((L, L, L, T))
    rng = Xoshiro(seed)
    g = hot_start(lat, rng)
    algs = (Heatbath(), Overrelaxation(), Overrelaxation(), Overrelaxation())
    println("# beta=$β L=$L T=$T rho=$ρ levels=$(join(levels, ",")) wlevel=$wlevel rmax=$rmax tmax=$tmax ",
            "therm=$ntherm skip=$nskip seed=$seed threads=$(Threads.nthreads())")
    for _ in 1:ntherm, alg in algs
        sweep!(g, β, alg, rng)
    end
    t0 = time()
    for c in 1:nconf
        for _ in 1:nskip, alg in algs
            sweep!(g, β, alg, rng)
        end
        O = zeros(length(levels), T)
        W = zeros(rmax, tmax)
        h, done = g, 0
        for n in all_levels
            h = stout_spatial(h, ρ, n - done)
            done = n
            k = findfirst(==(n), levels)
            k === nothing || (O[k, :] = slice_plaquettes(h))
            n == wlevel && (W = wilson_loops_rt(h, rmax, tmax))
        end
        buf = IOBuffer()
        @printf(buf, "%.12f", mean_plaquette(g))
        for k in 1:length(levels), t in 1:T
            @printf(buf, " %.12e", O[k, t])
        end
        for R in 1:rmax, τ in 1:tmax
            @printf(buf, " %.12e", W[R, τ])
        end
        println(String(take!(buf)))
        c % 50 == 0 && (flush(stdout); @printf(stderr, "beta=%.2f conf %d/%d, %.2f s per conf\n", β, c, nconf, (time() - t0) / c))
    end
end

main(ARGS)
