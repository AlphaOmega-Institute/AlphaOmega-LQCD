# 0⁺⁺ masses on the blocked fine lattice, for comparison with the coarse lattice in the same
# physical units: does the glueball have the same mass at the two lattice spacings?
#
#   julia --project=lattice lattice/examples/blocked_masses.jl fine.txt [--coarse coarse.txt] [options]
#
# fine.txt comes from examples/blocked_glueball.jl, coarse.txt from examples/glueball.jl.
# Options: --ops 4,8         stout levels (values) of the GEVP basis, on both lattices
#          --bin 10          configurations per bin
#          --t0 0            GEVP reference time, in each lattice's own units
#          --ratio 2.0015    a_coarse / a_fine (Necco–Sommer for β = 5.7 and 6.054)
# For each α and variant (fixed origin, zero momentum) of the fine file: the cosh effective
# mass of the GEVP ground-state correlator at fine times t = 0 … 7, in fine units and in coarse
# units (× ratio); the mass per fine step over the intervals 2j → 2j+2 (one mass per interval,
# from C(2j)/C(2j+2), so correlated adjacent ratios are not combined); and the squared overlap
# of each basis operator with the leading GEVP state of the window (not with the true glueball).
# With --coarse, the same for the coarse lattice at coarse times t = 0 … 3, and a table of
# matched intervals: the fine mass over 2j → 2j+2 against the coarse mass over j → j+1, per fine
# step, with their difference (the ensembles are independent, so the errors add in quadrature).
# Errors: delete-one jackknife over bins.

using LatticeYM, Printf, LinearAlgebra

fmt(x, e) = isfinite(x) && isfinite(e) ? @sprintf("%.3f(%d)", x, round(Int, e * 1e3)) : "   —    "

function operators(X, off, idx, T)
    O = zeros(size(X, 1), length(idx), T)
    for (j, k) in enumerate(idx), t in 1:T
        O[:, j, t] = X[:, off + (k - 1) * T + t]
    end
    return O
end

# effective masses of the GEVP ground state at t = 0 … tmax, the masses per step over the
# intervals lag·j → lag·(j+1) for j = 0 … jmax, then the squared overlaps with the GEVP state
function masses(O, b, t0, tmax, lag, jmax)
    P, M = correlator_products(O)
    Pb, Mb = bin_means(P, b), bin_means(M, b)
    n, T = size(O, 2), size(O, 3)
    function derived(ix)
        Pm = dropdims(sum(Pb[ix, :, :, :]; dims = 1); dims = 1) ./ length(ix)
        Mm = vec(sum(Mb[ix, :]; dims = 1)) ./ length(ix)
        C = correlator_from_means(Pm, Mm)
        _, V = gevp(C, t0, t0 + 1)
        c = [V[:, 1]' * C[:, :, t] * V[:, 1] for t in 1:T]
        f, _ = ground_state_fractions(C, collect(1:n); t0)
        mk = [cosh_mass(c[lag*j+1] / c[lag*(j+1)+1], lag * j, T; lag) for j in 0:jmax]
        return vcat(effective_masses(c)[1:tmax+1], mk, f)
    end
    return jackknife(derived, collect(1:size(Pb, 1)))
end

function report(est, err, tmax, jmax, lag, ops, scale, unit, other)
    row(label, s) = println(@sprintf("  %-30s", label), join([fmt(s * est[t+1], s * err[t+1]) for t in 0:tmax], " "))
    row("a m(t), $unit units, t = 0…$tmax", 1.0)
    row("the same in $other units", scale)
    if lag > 1
        println(@sprintf("  %-30s", "per step over $(lag)j → $(lag)j+$lag"),
                join([fmt(est[tmax+2+j], err[tmax+2+j]) for j in 0:jmax], " "))
    end
    o = tmax + jmax + 2
    println(@sprintf("  %-30s", "overlap² with GEVP state"),
            join(["level $(ops[j]) " * fmt(est[o+j], err[o+j]) for j in eachindex(ops)], "   "))
end

function main(args)
    fpath = args[1]
    opt = Dict("ops" => "4,8", "bin" => "10", "t0" => "0", "ratio" => "2.0015", "coarse" => "")
    i = 2
    while i <= length(args)
        key = lstrip(args[i], '-')
        haskey(opt, key) || error("unknown option $(args[i])")
        opt[key] = args[i+1]
        i += 2
    end
    ops = parse.(Int, split(opt["ops"], ","))
    b, t0 = parse(Int, opt["bin"]), parse(Int, opt["t0"])
    ratio = parse(Float64, opt["ratio"])

    _, kf, Xf, dropped = read_measurements(fpath)
    dropped && @warn "$fpath: ignoring an unfinished last line (file still being written?)"
    flev, Tf = parse.(Int, split(kf["levels"], ",")), parse(Int, kf["T"])
    alphas = parse.(Float64, split(kf["alphas"], ","))
    variants = haskey(kf, "variants") ? split(kf["variants"], ",") : ["fixed"]
    all(l -> l in flev, ops) || error("levels $ops must all be measured on the fine lattice ($flev)")
    size(Xf, 2) == 1 + length(alphas) * length(variants) * length(flev) * Tf ||
        error("fine rows have unexpected width $(size(Xf, 2))")
    idx = [findfirst(==(l), flev) for l in ops]
    tf = min(7, Tf ÷ 2 - 1)
    jf = min(3, (Tf ÷ 2) ÷ 2 - 1)
    matched = []
    println("fine: β = $(kf["beta"]), $(kf["L"])³ × $Tf blocked to $(kf["blocked"])³, $(size(Xf, 1)) configurations; ",
            "GEVP on stout levels $ops, t0 = $t0, bins of $b; a_c/a_f = $ratio")
    for (ia, α) in enumerate(alphas), (iv, var) in enumerate(variants)
        off = 1 + ((ia - 1) * length(variants) + iv - 1) * length(flev) * Tf
        est, err = masses(operators(Xf, off, idx, Tf), b, t0, tf, 2, jf)
        println("\nα = $α, $var")
        report(est, err, tf, jf, 2, ops, ratio, "fine", "coarse")
        push!(matched, ("α = $α, $var", est[tf+2:tf+2+jf], err[tf+2:tf+2+jf]))
    end

    isempty(opt["coarse"]) && return
    _, kc, Xc, dropc = read_measurements(opt["coarse"])
    dropc && @warn "$(opt["coarse"]): ignoring an unfinished last line"
    clev, Tc = parse.(Int, split(kc["levels"], ",")), parse(Int, kc["T"])
    all(l -> l in clev, ops) || error("levels $ops must all be measured on the coarse lattice ($clev)")
    tc = min(3, Tc ÷ 2 - 1)
    jc = min(3, Tc ÷ 2 - 1)
    est, err = masses(operators(Xc, 1, [findfirst(==(l), clev) for l in ops], Tc), b, t0, tc, 1, jc)
    println("\ncoarse: β = $(kc["beta"]), $(kc["L"])³ × $Tc, $(size(Xc, 1)) configurations")
    report(est, err, tc, jc, 1, ops, 1 / ratio, "coarse", "fine")

    J = min(jf, jc)
    mc, ec = est[tc+2:tc+2+J] ./ ratio, err[tc+2:tc+2+J] ./ ratio
    println("\nmatched intervals, mass per fine step (fine units): fine 2j → 2j+2 against coarse j → j+1")
    println(@sprintf("  %-26s", "j"), join([rpad("$j", 26) for j in 0:J]))
    println(@sprintf("  %-26s", "coarse"), join([rpad(fmt(mc[j+1], ec[j+1]), 26) for j in 0:J]))
    for (label, mf, ef) in matched
        cells = [rpad(fmt(mf[j+1], ef[j+1]) * " Δ " * fmt(mf[j+1] - mc[j+1], sqrt(ef[j+1]^2 + ec[j+1]^2)), 26) for j in 0:J]
        println(@sprintf("  %-26s", label), join(cells))
    end
end

main(ARGS)
