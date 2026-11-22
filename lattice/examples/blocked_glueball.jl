# Fine-lattice side of a comparison between two lattice spacings (Proposition 6 of
# docs/notes/spectral_transfer.tex). Each fine configuration is blocked once in the
# spatial directions with `block_spatial`, for several blocking parameters α and all
# eight blocking origins, and the 0⁺⁺ operators of examples/glueball.jl (stout levels
# with parameter ρ) are measured on each blocked lattice, on every fine time slice.
# Two variants are written for each α: the fixed-origin operator (origin 0), and the
# zero-momentum projection (operator amplitudes averaged over the eight origins).
# Correlators at fine lag 2t are compared with the coarse lattice at lag t.
#
#   julia --project=lattice -t auto lattice/examples/blocked_glueball.jl [options] > fine.txt
#
# Options: --beta 6.054      fine coupling (twice the r0/a of β = 5.7, Necco–Sommer)
#          --L 16 --T 32     fine extents; each blocked lattice is (L/2)³ × T
#          --therm 500 --confs 3000 --skip 4 --seed 6054
#          --alphas 0.0,0.3,0.6   blocking parameters
#          --rho 0.1 --levels 4,8,16,24   stout levels (a subset of the coarse ones)
#          --out file         write the output to this file instead of standard output
#          --checkpoint file  save the chain (fields and random state) every 25 configurations
#          --resume file      continue a chain from a checkpoint. With --out, the output file
#                             is first cut back to the configurations in the checkpoint, so
#                             lines written after it are not duplicated; then it is appended.
# Output: a header, then per configuration the fine mean plaquette followed, for each α,
# by the fixed-origin O[level, t] and then the zero-momentum O[level, t] (level-major).
# A resumed chain continues exactly as an uninterrupted one would.

using LatticeYM, Random, Printf, Serialization

function main(args)
    opt = Dict("beta" => "6.054", "L" => "16", "T" => "32", "therm" => "500", "confs" => "3000",
               "skip" => "4", "seed" => "6054", "alphas" => "0.0,0.3,0.6", "rho" => "0.1",
               "levels" => "4,8,16,24", "checkpoint" => "", "resume" => "", "out" => "")
    i = 1
    while i <= length(args)
        key = lstrip(args[i], '-')
        haskey(opt, key) || error("unknown option $(args[i])")
        opt[key] = args[i+1]
        i += 2
    end
    physics = Dict(k => opt[k] for k in ("beta", "L", "T", "therm", "skip", "seed", "alphas", "rho", "levels"))
    β, ρ = parse(Float64, opt["beta"]), parse(Float64, opt["rho"])
    L, T = parse(Int, opt["L"]), parse(Int, opt["T"])
    ntherm, nconf, nskip, seed = parse.(Int, (opt["therm"], opt["confs"], opt["skip"], opt["seed"]))
    alphas = parse.(Float64, split(opt["alphas"], ","))
    levels = parse.(Int, split(opt["levels"], ","))
    origins = vec([(a, b, c) for a in 0:1, b in 0:1, c in 0:1])
    algs = (Heatbath(), Overrelaxation(), Overrelaxation(), Overrelaxation())
    lat = Lattice((L, L, L, T))

    if isempty(opt["resume"])
        io = isempty(opt["out"]) ? stdout : open(opt["out"], "w")
        rng = Xoshiro(seed)
        g = hot_start(lat, rng)
        println(io, "# beta=$β L=$L T=$T blocked=$(L ÷ 2) rho=$ρ levels=$(join(levels, ",")) ",
                "alphas=$(join(alphas, ",")) variants=fixed,zero origins=8 therm=$ntherm skip=$nskip ",
                "seed=$seed threads=$(Threads.nthreads())")
        for _ in 1:ntherm, alg in algs
            sweep!(g, β, alg, rng)
        end
        done = 0
    else
        state = deserialize(opt["resume"])
        state.physics == physics || error("the checkpoint was made with different physics options: $(state.physics)")
        g, rng, done = GaugeField(lat, state.U), state.rng, state.done
        if isempty(opt["out"])
            io = stdout
            @printf(stderr, "resuming after configuration %d: keep only the header and the first %d data lines of the earlier output\n", done, done)
        else
            lines = readlines(opt["out"])
            length(lines) >= done + 1 || error("$(opt["out"]) has fewer lines than the checkpoint")
            open(opt["out"], "w") do f
                foreach(l -> println(f, l), lines[1:done+1])
            end
            io = open(opt["out"], "a")
            @printf(stderr, "resuming after configuration %d; %s cut back to it\n", done, opt["out"])
        end
    end

    t0 = time()
    for c in done+1:nconf
        for _ in 1:nskip, alg in algs
            sweep!(g, β, alg, rng)
        end
        buf = IOBuffer()
        @printf(buf, "%.12f", mean_plaquette(g))
        for α in alphas
            Ofix = zeros(length(levels), T)
            Ozero = zeros(length(levels), T)
            for o in origins
                O = glueball_operators(block_spatial(g, α; origin = o), ρ, levels)
                o == (0, 0, 0) && (Ofix .= O)
                Ozero .+= O ./ length(origins)
            end
            for O in (Ofix, Ozero), k in eachindex(levels), t in 1:T
                @printf(buf, " %.12e", O[k, t])
            end
        end
        println(io, String(take!(buf)))
        if c % 25 == 0
            flush(io)
            if !isempty(opt["checkpoint"])
                tmp = opt["checkpoint"] * ".tmp"
                serialize(tmp, (physics = physics, U = g.U, rng = copy(rng), done = c))
                mv(tmp, opt["checkpoint"]; force = true)
            end
            @printf(stderr, "beta=%.3f conf %d/%d, %.2f s per conf\n", β, c, nconf, (time() - t0) / (c - done))
        end
    end
    io === stdout || close(io)
end

main(ARGS)
