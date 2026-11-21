# Analyse the output of examples/glueball.jl: 0⁺⁺ effective masses per smearing
# level, a variational (GEVP) analysis, ground-state overlaps, and the string
# tension from smeared Wilson loops. Errors: jackknife over bins.
#
#   julia --project=lattice lattice/examples/glueball_analysis.jl data.txt [options]
#
# Options: --bin 10      configurations per bin
#          --drop 0      configurations dropped at the start
#          --ops 2,3,4   operator levels (1-based) used in the GEVP (default: all)
#          --t0 0        GEVP reference time (the eigenvectors are fixed at t0, t0+1)
#          --tau 2       Wilson-loop time used for V(R) = log W(R,τ)/W(R,τ+1)

using LatticeYM, Printf, LinearAlgebra

function read_data(path)
    header = ""
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        isempty(strip(line)) && continue
        if startswith(line, "#")
            header = line
        else
            push!(rows, parse.(Float64, split(line)))
        end
    end
    kv = Dict(split(p, "=")[1] => split(p, "=")[2] for p in split(header[3:end]) if occursin("=", p))
    return kv, permutedims(reduce(hcat, rows))
end

fmt(x, e) = isfinite(x) && isfinite(e) ? @sprintf("%.4f(%s)", x, lpad(round(Int, e * 1e4), 1)) : "   —   "

function main(args)
    path = args[1]
    opt = Dict("bin" => "10", "drop" => "0", "ops" => "", "t0" => "0", "tau" => "2")
    i = 2
    while i <= length(args)
        opt[lstrip(args[i], '-')] = args[i+1]
        i += 2
    end
    kv, X = read_data(path)
    T, β = parse(Int, kv["T"]), parse(Float64, kv["beta"])
    levels = parse.(Int, split(kv["levels"], ","))
    rmax, tmax = parse(Int, kv["rmax"]), parse(Int, kv["tmax"])
    nlev = length(levels)
    b, drop, t0, τ = parse.(Int, (opt["bin"], opt["drop"], opt["t0"], opt["tau"]))
    ops = isempty(opt["ops"]) ? collect(1:nlev) : parse.(Int, split(opt["ops"], ","))
    X = X[drop+1:end, :]
    nconf = size(X, 1)

    plaq = X[:, 1]
    O = zeros(nconf, nlev, T)
    for k in 1:nlev, t in 1:T
        O[:, k, t] = X[:, 1 + (k - 1) * T + t]
    end
    W = zeros(nconf, rmax, tmax)
    for R in 1:rmax, s in 1:tmax
        W[:, R, s] = X[:, 1 + nlev * T + (R - 1) * tmax + s]
    end

    P, M = correlator_products(O)
    Pb, Mb, Wb = bin_means(P, b), bin_means(M, b), bin_means(W, b)
    nb = size(Pb, 1)
    tshow = min(4, T ÷ 2 - 1)

    # everything derived from one set of bins, as a flat vector
    function derived(ix)
        Pm = dropdims(sum(Pb[ix, :, :, :]; dims = 1); dims = 1) ./ length(ix)
        Mm = vec(sum(Mb[ix, :]; dims = 1)) ./ length(ix)
        Wm = dropdims(sum(Wb[ix, :, :]; dims = 1); dims = 1) ./ length(ix)
        C = correlator_from_means(Pm, Mm)
        out = Float64[]
        for k in 1:nlev                                     # effective masses per level
            append!(out, effective_masses(C[k, k, :])[1:tshow+1])
        end
        Cs = C[ops, ops, :]
        λ, V = gevp(Cs, t0, t0 + 1)
        for n in 1:min(2, length(ops))                      # projected correlators
            c = [V[:, n]' * Cs[:, :, t] * V[:, n] for t in 1:T]
            append!(out, effective_masses(c)[1:tshow+1])
        end
        length(ops) == 1 && append!(out, fill(NaN, tshow + 1))
        # ground-state overlap of each level: (C(t0) v0)_k^2 / C_kk(0) e^{-2 E0 t0}
        v0 = zeros(nlev); v0[ops] = V[:, 1]
        c0 = [v0' * C[:, :, t] * v0 for t in 1:T]
        E0 = cosh_mass(c0[t0+2] / c0[t0+3], t0 + 1, T)
        Zt0 = C[:, :, t0+1] * v0 ./ sqrt(v0' * C[:, :, t0+1] * v0)
        f = [Zt0[k]^2 * exp(2E0 * t0) / C[k, k, 1] for k in 1:nlev]
        append!(out, f)
        # first excited level in the channel (GEVP n = 1), and the defect bound of
        # Proposition 1(ii) for the one-dimensional window of each level:
        # δ ≥ (Λ − E0) √(1/f − 1), with Λ ≈ E1
        E1 = NaN
        if length(ops) >= 2
            v1 = zeros(nlev); v1[ops] = V[:, 2]
            c1 = [v1' * C[:, :, t] * v1 for t in 1:T]
            E1 = cosh_mass(c1[t0+1] / c1[t0+2], t0, T)
        end
        append!(out, [(E1 - E0) * sqrt(max(1 / fk - 1, 0.0)) for fk in f])
        # second estimator, from the decay of each diagonal correlator at t = 1 and 2:
        # f_k(t) = C_kk(t) cosh(E0 T/2) / (C_kk(0) cosh(E0 (t − T/2))); excited states make it
        # an overestimate that decreases with t
        for t in (1, 2), k in 1:nlev
            push!(out, C[k, k, t+1] * cosh(E0 * T / 2) / (C[k, k, 1] * cosh(E0 * (t - T / 2))))
        end
        # static potential and Cornell fits
        Vr = [log(Wm[R, τ] / Wm[R, τ+1]) for R in 1:rmax]
        append!(out, Vr)
        Rs = collect(2.0:rmax)
        A3 = hcat(ones(length(Rs)), Rs, -1 ./ Rs)
        p3 = A3 \ Vr[2:end]                                 # V0, σ, e free
        A2 = hcat(ones(length(Rs)), Rs)
        p2 = A2 \ (Vr[2:end] .+ (π / 12) ./ Rs)             # e = π/12 fixed
        σ = p2[2]
        append!(out, [p3[2], p3[3], σ, E0, E0 / sqrt(abs(σ))])
        return out
    end

    est, err = jackknife(derived, collect(1:nb))
    pm, pe, pτ = mean_err(plaq)
    @printf("β = %.3f, %s, %d configurations in %d bins of %d, ⟨P⟩ = %.6f(%d), τ_int = %.1f\n",
            β, "L=$(kv["L"]) T=$T", nconf, nb, b, pm, round(Int, pe * 1e6), pτ)
    println("smearing levels: ", levels, "  (ρ = ", kv["rho"], ")   GEVP operators: ", levels[ops], "  t0 = $t0")
    p = 0
    println("\neffective mass a·m_eff(t) per smearing level, t = 0 … $tshow")
    for k in 1:nlev
        println(@sprintf("  level %3d: ", levels[k]), join([fmt(est[p+j], err[p+j]) for j in 1:tshow+1], "  "))
        p += tshow + 1
    end
    for n in 1:2
        println(@sprintf("  GEVP n=%d : ", n - 1), join([fmt(est[p+j], err[p+j]) for j in 1:tshow+1], "  "))
        p += tshow + 1
    end
    println("\nground-state fraction f_k = |<G|O_k|0>|²/C_kk(0) = ‖Pφ₀‖² per level")
    println("  ", join([@sprintf("%d: %s", levels[k], fmt(est[p+k], err[p+k])) for k in 1:nlev], "   "))
    p += nlev
    println("defect bound a·δ_k ≥ a(E₁ − E₀)·√(1/f_k − 1) per level (Proposition 1(ii), Λ ≈ E₁)")
    println("  ", join([@sprintf("%d: %s", levels[k], fmt(est[p+k], err[p+k])) for k in 1:nlev], "   "))
    p += nlev
    for t in (1, 2)
        println("ground-state fraction from the decay at t = $t (overestimate)")
        println("  ", join([@sprintf("%d: %s", levels[k], fmt(est[p+k], err[p+k])) for k in 1:nlev], "   "))
        p += nlev
    end
    println("\nstatic potential a·V(R) from τ = $τ → $(τ+1), R = 1 … $rmax")
    println("  ", join([fmt(est[p+R], err[p+R]) for R in 1:rmax], "  "))
    p += rmax
    labels = ("σa² (Cornell, e free)", "e (Cornell, e free)", "σa² (e = π/12)", "a·m₀ (GEVP, t0+1)", "m₀/√σ")
    println()
    for (j, l) in enumerate(labels)
        println(@sprintf("  %-22s %s", l, fmt(est[p+j], err[p+j])))
    end
end

main(ARGS)
