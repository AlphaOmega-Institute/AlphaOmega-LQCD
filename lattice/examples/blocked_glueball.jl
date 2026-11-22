# Fine-lattice side of a comparison between two lattice spacings (Proposition 6 of
# docs/notes/spectral_transfer.tex). Each fine configuration is blocked once in the
# spatial directions with `block_spatial`, for several blocking parameters α and all
# eight blocking origins, and the 0⁺⁺ operators of examples/glueball.jl (stout levels
# with parameter ρ) are measured on each blocked lattice, on every fine time slice.
# Two variants are written for each α: the fixed-origin operator (origin 0), and the
# zero-momentum projection (operator amplitudes averaged over the eight origins).
# Correlators at fine lag 2t are compared with the coarse lattice at lag t.
#
#   julia --project=lattice -t auto lattice/examples/blocked_glueball.jl --out fine.txt --checkpoint fine.ckpt
#
# Physics options (recorded in the header and the checkpoint):
#          --beta 6.054      fine coupling (twice the r0/a of β = 5.7, Necco–Sommer)
#          --L 16 --T 32     fine extents; each blocked lattice is (L/2)³ × T
#          --therm 500 --skip 4 --seed 6054
#          --alphas 0.0,0.3,0.6   blocking parameters
#          --rho 0.1 --levels 4,8,16,24   stout levels (a subset of the coarse ones)
# Run options:
#          --confs 3000      total number of configurations (on resume: the new total)
#          --out file        output file; without it the output goes to standard output
#                            and checkpoints are not available
#          --checkpoint file save the chain every 25 configurations and at the last one:
#                            fields, random state, and the length and FNV-1a hash of the
#                            output written so far
#          --resume file     continue the chain of a checkpoint (checkpoints then go on
#                            to the same file unless --checkpoint names another). Physics
#                            options come from the checkpoint; any given explicitly must
#                            agree with it. The output file must begin with exactly the
#                            bytes the checkpoint recorded (same header, same rows);
#                            anything after them, written after the checkpoint, is cut
#                            off in place and the run appends from there.
# A fresh run refuses to overwrite an existing output or checkpoint file.
# Output: a header, then per configuration the fine mean plaquette followed, for each α,
# by the fixed-origin O[level, t] and then the zero-momentum O[level, t] (level-major).
# A resumed chain continues exactly as an uninterrupted one would, with any number of
# threads (the updates draw a fixed number of per-chunk seeds).

using LatticeYM, Random, Printf, Serialization

const PHYSICS = ("beta", "L", "T", "therm", "skip", "seed", "alphas", "rho", "levels")

function sync!(io)
    flush(io)
    Sys.isunix() && ccall(:fsync, Cint, (Cint,), fd(io))
    return nothing
end

function save_checkpoint(path, state)
    tmp = path * ".tmp"
    open(tmp, "w") do f
        serialize(f, state)
        sync!(f)
    end
    mv(tmp, path; force = true)
end

function main(args)
    opt = Dict("beta" => "6.054", "L" => "16", "T" => "32", "therm" => "500", "confs" => "3000",
               "skip" => "4", "seed" => "6054", "alphas" => "0.0,0.3,0.6", "rho" => "0.1",
               "levels" => "4,8,16,24", "checkpoint" => "", "resume" => "", "out" => "")
    given = Set{String}()
    i = 1
    while i <= length(args)
        key = lstrip(args[i], '-')
        haskey(opt, key) || error("unknown option $(args[i])")
        i < length(args) || error("option $(args[i]) needs a value")
        opt[key] = args[i+1]
        push!(given, key)
        i += 2
    end
    out, ckpt, resume = opt["out"], opt["checkpoint"], opt["resume"]
    isempty(ckpt) && isempty(resume) || !isempty(out) ||
        error("--checkpoint and --resume need --out: the checkpoint records the committed part of that file")

    state = nothing
    if !isempty(resume)
        state = deserialize(resume)
        hasproperty(state, :nbytes) ||
            error("$resume was written by an older version of this script, which did not record the output; start a new run")
        for k in PHYSICS
            k in given && opt[k] != state.physics[k] &&
                error("--$k $(opt[k]) differs from the checkpoint ($(state.physics[k]))")
            opt[k] = state.physics[k]
        end
        isempty(ckpt) && (ckpt = resume)
    else
        isempty(out) || !ispath(out) || error("$out exists; refusing to overwrite it (resume it with its checkpoint, or choose another name)")
        isempty(ckpt) || !ispath(ckpt) || error("$ckpt exists; refusing to overwrite it")
    end
    physics = Dict(k => opt[k] for k in PHYSICS)
    β, ρ = parse(Float64, opt["beta"]), parse(Float64, opt["rho"])
    L, T = parse(Int, opt["L"]), parse(Int, opt["T"])
    ntherm, nconf, nskip, seed = parse.(Int, (opt["therm"], opt["confs"], opt["skip"], opt["seed"]))
    alphas = parse.(Float64, split(opt["alphas"], ","))
    levels = parse.(Int, split(opt["levels"], ","))
    origins = vec([(a, b, c) for a in 0:1, b in 0:1, c in 0:1])
    algs = (Heatbath(), Overrelaxation(), Overrelaxation(), Overrelaxation())
    lat = Lattice((L, L, L, T))

    if state === nothing
        io = isempty(out) ? stdout : open(out, "w")
        rng = Xoshiro(seed)
        g = hot_start(lat, rng)
        header = "# beta=$β L=$L T=$T blocked=$(L ÷ 2) rho=$ρ levels=$(join(levels, ",")) " *
                 "alphas=$(join(alphas, ",")) variants=fixed,zero origins=8 therm=$ntherm skip=$nskip " *
                 "seed=$seed threads=$(Threads.nthreads())"
        write(io, header, "\n")
        hash, nbytes = fnv1a(header * "\n"), sizeof(header) + 1
        for _ in 1:ntherm, alg in algs
            sweep!(g, β, alg, rng)
        end
        done = 0
    else
        header, hash, nbytes, done = state.header, state.hash, state.nbytes, state.done
        isfile(out) || error("$out not found")
        size0 = filesize(out)
        size0 >= nbytes || error("$out has $size0 bytes, fewer than the $nbytes recorded at configuration $done")
        firstline = open(readline, out)
        firstline == header || error("$out has a different header from the checkpoint (another run?):\n  file:       $firstline\n  checkpoint: $header")
        file_prefix_hash(out, nbytes) == hash ||
            error("the first $nbytes bytes of $out differ from those recorded in the checkpoint (another run, or a modified file)")
        open(out, "r+") do f
            truncate(f, nbytes)
        end
        io = open(out, "a")
        g, rng = GaugeField(lat, state.U), state.rng
        LatticeYM.SU3_FALLBACKS[] = state.fallbacks
        @printf(stderr, "resuming after configuration %d: %s verified, cut from %d to %d bytes\n", done, out, size0, nbytes)
        done >= nconf && @printf(stderr, "the checkpoint already holds %d ≥ %d configurations; nothing to do\n", done, nconf)
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
        write(buf, '\n')
        line = take!(buf)
        write(io, line)
        hash, nbytes = fnv1a(line, hash), nbytes + length(line)
        if c % 25 == 0 || c == nconf
            flush(io)
            if !isempty(ckpt)
                sync!(io)
                filesize(out) == nbytes ||
                    error("$out has $(filesize(out)) bytes where $nbytes were written: is another process writing to it?")
                save_checkpoint(ckpt, (physics = physics, U = g.U, rng = copy(rng), done = c, header = header,
                                       nbytes = nbytes, hash = hash, fallbacks = LatticeYM.SU3_FALLBACKS[]))
            end
            @printf(stderr, "beta=%.3f conf %d/%d, %.2f s per conf, projection fallbacks %d\n",
                    β, c, nconf, (time() - t0) / (c - done), LatticeYM.SU3_FALLBACKS[])
        end
    end
    io === stdout || close(io)
end

main(ARGS)
