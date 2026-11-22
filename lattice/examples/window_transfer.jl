# Compare a coarse lattice with a blocked fine lattice through an operator window
# (Proposition 6 of docs/notes/spectral_transfer.tex).
#
#   julia --project=lattice lattice/examples/window_transfer.jl coarse.txt fine.txt [options]
#
# coarse.txt comes from examples/glueball.jl, fine.txt from examples/blocked_glueball.jl.
# The fine file must use a subset of the coarse smearing levels and the same ρ, a time
# extent twice the coarse one, and blocked spatial extent equal to the coarse one.
# Options:
#          --window 4,8      stout levels spanning the window (values, not indices)
#          --bin 10          configurations per bin, on both lattices
#          --drop 0          fine configurations dropped at the start
# A file still being written can be analysed: only lines ending in a newline are read
# (`read_measurements`), and every row read is checked.
# For each blocking parameter α and each variant (fixed origin, zero momentum) the script
# reports the coarse compression K_c and leakage L_c at coarse lag 1, the fine ones K_f and
# L_f at fine lag 2, the alignment V induced by the correspondence O ↦ O∘B (polar part of
# G_f^{1/2} G_c^{-1/2}), the normalization mismatch ‖M − I‖, the in-window mismatch
# ‖K_f V − V K_c‖, the larger leakage, the bound of Proposition 6, and the defect itself
# (`window_defect`). It then reads the window map in the doubled space H_c ⊕ H_f
# (`graph_departure`): the Krein square ‖M² − I‖ of the raw correspondence's graph, and
# how far one step of evolution carries the aligned graph of V into the anti-aligned graph
# (half the in-window mismatch) and out of both windows, and in total. The operators are
# centered (vacuum expectation subtracted) on each lattice, so the map compared is the one
# between centered excitations. Errors: delete-one jackknife over the
# bins of each lattice, the two independent variances added; all α and both variants are
# computed from the same fine bins, so they are resampled together. When both variants
# are present, the differences fixed − zero of the mismatch, the bound and the defect are
# reported for each α with paired errors (same configurations, same bins).

using LatticeYM, Printf, LinearAlgebra

function read_file(path)
    _, kv, X, dropped = read_measurements(path)
    dropped && @warn "$path: ignoring an unfinished last line (file still being written?)"
    return kv, X
end

binned_products(O, b) = ((P, M) = correlator_products(O); (bin_means(P, b), bin_means(M, b)))
meanC(Pb, Mb, ix) = correlator_from_means(dropdims(sum(Pb[ix, :, :, :]; dims = 1); dims = 1) ./ length(ix),
                                          vec(sum(Mb[ix, :]; dims = 1)) ./ length(ix))

function window_quantities(Cc, Cf)
    Kc, Lc, Gc = window_compression(Cc, 1)
    Kf, Lf, Gf = window_compression(Cf, 2)
    V, M = induced_alignment(Gf, Gc)
    λc, λf = eigvals(Symmetric(Lc)), eigvals(Symmetric(Lf))
    leak = sqrt(max(maximum(λc), maximum(λf), 0.0))
    mis = opnorm(Kf * V - V * Kc)
    anti, out, tot = graph_departure(Kf, Lf, Kc, Lc, V)
    return vcat(sort(eigvals(Symmetric(Kc)); rev = true), sort(eigvals(Symmetric(Kf)); rev = true),
                sort(λc), sort(λf),
                [opnorm(M - I), mis, leak, mis + leak, window_defect(Kf, Lf, Kc, Lc, V),
                 opnorm(M^2 - I), anti, out, tot])
end

function two_sample_jackknife(f, nc, nf)
    full = f(1:nc, 1:nf)
    var = zero(full)
    for (n, g) in ((nc, i -> f(setdiff(1:nc, i), 1:nf)), (nf, i -> f(1:nc, setdiff(1:nf, i))))
        s = [g(i) for i in 1:n]
        m = sum(s) / n
        var = var .+ sum(x -> (x .- m) .^ 2, s) .* ((n - 1) / n)
    end
    return full, sqrt.(var)
end

fmt(x, e) = @sprintf("%.4f(%d)", x, round(Int, e * 1e4))

function main(args)
    cpath, fpath = args[1], args[2]
    opt = Dict("window" => "4,8", "bin" => "10", "drop" => "0")
    i = 3
    while i <= length(args)
        key = lstrip(args[i], '-')
        haskey(opt, key) || error("unknown option $(args[i])")
        opt[key] = args[i+1]
        i += 2
    end
    wlev = parse.(Int, split(opt["window"], ","))
    b, drop = parse(Int, opt["bin"]), parse(Int, opt["drop"])
    kc, Xc = read_file(cpath)
    kf, Xf = read_file(fpath)
    Xf = Xf[drop+1:end, :]

    clev, flev = parse.(Int, split(kc["levels"], ",")), parse.(Int, split(kf["levels"], ","))
    Lc, Tc, Lf, Tf = parse.(Int, (kc["L"], kc["T"], kf["L"], kf["T"]))
    rmax, tmax = parse(Int, kc["rmax"]), parse(Int, kc["tmax"])
    alphas = parse.(Float64, split(kf["alphas"], ","))
    variants = haskey(kf, "variants") ? split(kf["variants"], ",") : ["fixed"]
    kc["rho"] == kf["rho"] || error("different stout parameters: ρ = $(kc["rho"]) and $(kf["rho"])")
    Lf == 2Lc || error("the blocked fine lattice ($(Lf ÷ 2)³) does not match the coarse one ($(Lc)³)")
    Tf == 2Tc || error("the fine time extent ($Tf) is not twice the coarse one ($Tc)")
    all(l -> l in clev && l in flev, wlev) || error("window levels $wlev must be measured on both lattices")
    size(Xc, 2) == 1 + length(clev) * Tc + rmax * tmax || error("coarse rows have unexpected width $(size(Xc, 2))")
    nvar = length(alphas) * length(variants)
    size(Xf, 2) == 1 + nvar * length(flev) * Tf || error("fine rows have unexpected width $(size(Xf, 2))")

    wc = [findfirst(==(l), clev) for l in wlev]
    wf = [findfirst(==(l), flev) for l in wlev]
    Oc = zeros(size(Xc, 1), length(wc), Tc)
    for (j, k) in enumerate(wc), t in 1:Tc
        Oc[:, j, t] = Xc[:, 1 + (k - 1) * Tc + t]
    end
    Pc, Mc = binned_products(Oc, b)
    fine = map(1:nvar) do v
        Of = zeros(size(Xf, 1), length(wf), Tf)
        off = 1 + (v - 1) * length(flev) * Tf
        for (j, k) in enumerate(wf), t in 1:Tf
            Of[:, j, t] = Xf[:, off + (k - 1) * Tf + t]
        end
        binned_products(Of, b)
    end
    nfb = size(fine[1][1], 1)

    n = length(wlev)
    per = 4n + 9
    paired = variants == ["fixed", "zero"]
    pix = [4n+2, 4n+4, 4n+5, 4n+9]  # mismatch, bound, defect, total departure from the aligned graph
    function f(ic, jf)
        q = [window_quantities(meanC(Pc, Mc, ic), meanC(Pf, Mf, jf)) for (Pf, Mf) in fine]
        v = vcat(q...)
        if paired                     # fixed − zero
            for ia in eachindex(alphas)
                append!(v, q[2ia-1][pix] .- q[2ia][pix])
            end
        end
        return v
    end
    est, err = two_sample_jackknife(f, size(Pc, 1), nfb)

    println("coarse: β = $(kc["beta"]), $(Lc)³ × $Tc, $(size(Xc, 1)) configurations; ",
            "fine: β = $(kf["beta"]), $(Lf)³ × $Tf blocked to $(Lf ÷ 2)³, $(size(Xf, 1)) configurations")
    println("window: stout levels $wlev (ρ = $(kc["rho"])), bins of $b; centered excitations")
    for (ia, α) in enumerate(alphas), (iv, var) in enumerate(variants)
        p = ((ia - 1) * length(variants) + (iv - 1)) * per
        println("\nα = $α, $var")
        for (label, m) in (("eig K_c (coarse lag 1)", n), ("eig K_f (fine lag 2)", n), ("eig L_c", n), ("eig L_f", n))
            println(@sprintf("  %-24s", label), join([fmt(est[p+j], err[p+j]) for j in 1:m], "  "))
            p += m
        end
        for label in ("‖M − I‖ (normalization)", "‖K_f V − V K_c‖", "max leakage^{1/2}",
                      "bound (Proposition 6)", "defect (window_defect)",
                      "Krein square ‖M² − I‖", "graph: to anti-aligned", "graph: out of windows",
                      "graph: total departure")
            p += 1
            println(@sprintf("  %-24s %s", label, fmt(est[p], err[p])))
        end
    end
    if paired
        println("\nfixed origin − zero momentum (paired errors: same configurations and bins)")
        for (ia, α) in enumerate(alphas)
            p = nvar * per + length(pix) * (ia - 1)
            println("  α = $α   Δ‖K_f V − V K_c‖ ", fmt(est[p+1], err[p+1]), "   Δbound ", fmt(est[p+2], err[p+2]),
                    "   Δdefect ", fmt(est[p+3], err[p+3]), "   Δdeparture ", fmt(est[p+4], err[p+4]))
        end
    end
end

main(ARGS)
