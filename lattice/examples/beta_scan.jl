# Average plaquette and action density of SU(3) Yang–Mills theory at a list of couplings.
#
#   julia --project=lattice -t auto lattice/examples/beta_scan.jl [options] β₁ β₂ ...
#
# Options: --L 8        lattice extent (L⁴)
#          --twist      symmetric 't Hooft twist in the (1,2) and (3,4) planes
#          --therm 200  thermalization sweeps
#          --sweeps 2000 measured sweeps
#          --seed 1
# One sweep is one heatbath pass followed by three overrelaxation passes.

using LatticeYM, Random, Printf

function main(args)
    L, twist, ntherm, nmeas, seed = 8, (), 200, 2000, 1
    βs = Float64[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--twist"
            twist = ((1, 2) => 1, (3, 4) => 1)
            i += 1
        elseif a in ("--L", "--therm", "--sweeps", "--seed")
            v = parse(Int, args[i+1])
            a == "--L" ? (L = v) : a == "--therm" ? (ntherm = v) : a == "--sweeps" ? (nmeas = v) : (seed = v)
            i += 2
        else
            push!(βs, parse(Float64, a))
            i += 1
        end
    end
    isempty(βs) && error("give at least one β")

    lat = Lattice((L, L, L, L); twist)
    algs = (Heatbath(), Overrelaxation(), Overrelaxation(), Overrelaxation())
    println("# SU(3) Wilson action on $(L)^4, twist = ", isempty(twist) ? "none" : twist,
            ", $ntherm + $nmeas sweeps, $(Threads.nthreads()) threads")
    println("#     β        ⟨P⟩                 ⟨S⟩ = 1 − ⟨P⟩         τ_int")
    for β in βs
        rng = Xoshiro(seed)
        g = hot_start(lat, rng)
        for _ in 1:ntherm, alg in algs
            sweep!(g, β, alg, rng)
        end
        P = map(1:nmeas) do _
            foreach(alg -> sweep!(g, β, alg, rng), algs)
            mean_plaquette(g)
        end
        m, e, τ = mean_err(P)
        @printf("%8.3f   %.6f ± %.6f   %.6f ± %.6f   %5.2f\n", β, m, e, 1 - m, e, τ)
    end
end

main(ARGS)
