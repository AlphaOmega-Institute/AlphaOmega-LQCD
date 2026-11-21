# Fine-lattice side of a comparison between two lattice spacings (Proposition 6 of
# docs/notes/spectral_transfer.tex). Each fine configuration is blocked once in the
# spatial directions with `block_spatial`, for several blocking parameters α, and
# the same 0⁺⁺ operators as in examples/glueball.jl (stout levels with parameter ρ)
# are measured on the blocked lattice, on every fine time slice. Correlators at fine
# lag 2t are then compared with the coarse lattice at lag t.
#
#   julia --project=lattice -t auto lattice/examples/blocked_glueball.jl [options] > fine.txt
#
# Options: --beta 6.054      fine coupling (twice the r0/a of β = 5.7, Necco–Sommer)
#          --L 16 --T 32     fine extents; the blocked lattice is (L/2)³ × T
#          --therm 500 --confs 1000 --skip 4 --seed 1
#          --alphas 0.0,0.3,0.6   blocking parameters
#          --rho 0.1 --levels 0,4,8,16,24,32   as for the coarse lattice
# Output: a header, then per configuration the fine mean plaquette followed, for each
# α in turn, by O[level, t] (level-major) on the blocked lattice.

using LatticeYM, Random, Printf

function main(args)
    opt = Dict("beta" => "6.054", "L" => "16", "T" => "32", "therm" => "500", "confs" => "1000",
               "skip" => "4", "seed" => "1", "alphas" => "0.0,0.3,0.6", "rho" => "0.1",
               "levels" => "0,4,8,16,24,32")
    i = 1
    while i <= length(args)
        key = lstrip(args[i], '-')
        haskey(opt, key) || error("unknown option $(args[i])")
        opt[key] = args[i+1]
        i += 2
    end
    β, ρ = parse(Float64, opt["beta"]), parse(Float64, opt["rho"])
    L, T = parse(Int, opt["L"]), parse(Int, opt["T"])
    ntherm, nconf, nskip, seed = parse.(Int, (opt["therm"], opt["confs"], opt["skip"], opt["seed"]))
    alphas = parse.(Float64, split(opt["alphas"], ","))
    levels = parse.(Int, split(opt["levels"], ","))

    rng = Xoshiro(seed)
    g = hot_start(Lattice((L, L, L, T)), rng)
    algs = (Heatbath(), Overrelaxation(), Overrelaxation(), Overrelaxation())
    println("# beta=$β L=$L T=$T blocked=$(L ÷ 2) rho=$ρ levels=$(join(levels, ",")) ",
            "alphas=$(join(alphas, ",")) therm=$ntherm skip=$nskip seed=$seed threads=$(Threads.nthreads())")
    for _ in 1:ntherm, alg in algs
        sweep!(g, β, alg, rng)
    end
    t0 = time()
    for c in 1:nconf
        for _ in 1:nskip, alg in algs
            sweep!(g, β, alg, rng)
        end
        buf = IOBuffer()
        @printf(buf, "%.12f", mean_plaquette(g))
        for α in alphas
            O = glueball_operators(block_spatial(g, α), ρ, levels)
            for k in eachindex(levels), t in 1:T
                @printf(buf, " %.12e", O[k, t])
            end
        end
        println(String(take!(buf)))
        c % 25 == 0 && (flush(stdout); @printf(stderr, "beta=%.3f conf %d/%d, %.2f s per conf\n", β, c, nconf, (time() - t0) / c))
    end
end

main(ARGS)
