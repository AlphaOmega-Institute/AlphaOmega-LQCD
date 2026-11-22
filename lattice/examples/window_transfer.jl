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
#          --lag 1           coarse lag ℓ; the fine lag is 2ℓ (ℓ = 2 compares over twice the
#                            physical time; the compressions and leakages are recomputed from the
#                            correlators at those lags, K(2τ) = K(τ)² + L(τ), not by squaring K)
#          --mref m          reference mass for the one-sided bound, per coarse step, in coarse
#                            units: s = e^{−ℓ m} (default: the top coarse window eigenvalue,
#                            s = λ_max(K_c)). With --mref, the defect and ε₊ are also printed
#                            relative to the reference gap 1 − s
#          --qgrid 0.1,0.2,0.3   omitted tops q at which the one-sided bound μ(q) is printed
# A file still being written can be analysed: only lines ending in a newline are read
# (`read_measurements`), and every row read is checked.
# For each blocking parameter α and each variant (fixed origin, zero momentum) the script
# reports the coarse compression K_c and leakage L_c at coarse lag ℓ (default 1), the fine ones K_f and
# L_f at fine lag 2ℓ, the alignment V induced by the correspondence O ↦ O∘B (polar part of
# G_f^{1/2} G_c^{-1/2}), the normalization mismatch ‖M − I‖, the in-window mismatch
# ‖K_f V − V K_c‖, the larger leakage, the bound of Proposition 6, and the defect itself
# (`window_defect`). It then reads the window map in the doubled space H_c ⊕ H_f
# (`graph_departure`): ‖M² − I‖, the largest absolute Krein square (squared-length imbalance)
# of the raw correspondence's graph, and
# how far one step of evolution carries the aligned graph of V into the anti-aligned graph
# (half the in-window mismatch) and out of both windows, and in total. The operators are
# centered (vacuum expectation subtracted) on each lattice, so the map compared is the one
# between centered excitations.
# One-sided diagnostic (docs/notes/README.md, "One-sided bound"): the signed extremal eigenvalues
# of V*K_f V − K_c and ε₊ = max(0, λ_max), the only part that can make the fine evolution decay
# more slowly than the coarse one; the floor on the omitted top q = ‖QTQ‖ that positivity of the
# fine transfer matrix implies, q ≥ λ_max(K_f^{-1/2} L_f K_f^{-1/2}); the largest q for which the
# bound μ = top of [k ℓ; ℓ q] stays below the reference s, with k the fine window top measured
# directly or inferred from the coarse lattice as λ_max(K_c) + ε₊ (ℓ² = λ_max(L_f)); and μ(q) on
# the grid for both choices of k. Keeping the direction of the one-sided error sharpens the
# coarse-derived bound: V*K_f V = K_c + D ≤ K_c + D₊, so k ≤ λ_max(K_c + D₊) ≤ λ_max(K_c) + ε₊.
# Also printed: the angle between the leading eigenvectors of K_c and of V*K_f V, the rotation
# that the blocking-fixed pairing introduces relative to a pairing in spectral order. All are
# computed in every jackknife replicate. Errors: delete-one jackknife over the
# bins of each lattice, the two independent variances added; all α and both variants are
# computed from the same fine bins, so they are resampled together. When both variants
# are present, the differences fixed − zero of the mismatch, the bound, the defect, the total
# departure and ε₊ are reported for each α with paired errors (same configurations, same bins).

using LatticeYM, Printf, LinearAlgebra

function read_file(path)
    _, kv, X, dropped = read_measurements(path)
    dropped && @warn "$path: ignoring an unfinished last line (file still being written?)"
    return kv, X
end

binned_products(O, b) = ((P, M) = correlator_products(O); (bin_means(P, b), bin_means(M, b)))
meanC(Pb, Mb, ix) = correlator_from_means(dropdims(sum(Pb[ix, :, :, :]; dims = 1); dims = 1) ./ length(ix),
                                          vec(sum(Mb[ix, :]; dims = 1)) ./ length(ix))

function window_quantities(Cc, Cf, sref, qgrid, lag)
    Kc, Lc, Gc = window_compression(Cc, lag)
    Kf, Lf, Gf = window_compression(Cf, 2lag)
    V, M = induced_alignment(Gf, Gc)
    λc, λf = eigvals(Symmetric(Lc)), eigvals(Symmetric(Lf))
    leak = sqrt(max(maximum(λc), maximum(λf), 0.0))
    mis = opnorm(Kf * V - V * Kc)
    anti, out, tot = graph_departure(Kf, Lf, Kc, Lc, V)
    dl = eigvals(Hermitian((V' * Kf * V - Kc + (V' * Kf * V - Kc)') / 2))
    εp = max(dl[end], 0.0)
    kd, kc = eigmax(Hermitian(Kf)), eigmax(Hermitian(Kc)) + εp
    ℓ2 = max(maximum(λf), 0.0)
    s = sref === nothing ? eigmax(Hermitian(Kc)) : sref
    FD = eigen(Hermitian((V' * Kf * V - Kc + (V' * Kf * V - Kc)') / 2))
    Dp = FD.vectors * Diagonal(max.(FD.values, 0.0)) * FD.vectors'
    km = eigmax(Hermitian(Kc + Dp))
    θ = acosd(min(1.0, abs(dot(eigen(Hermitian(Kc)).vectors[:, end], V' * eigen(Hermitian(Kf)).vectors[:, end]))))
    def = window_defect(Kf, Lf, Kc, Lc, V)
    rel = sref === nothing ? [NaN, NaN] : [def, εp] ./ (1 - s)
    qmd = q_needed_matrix(s, Kf, Lf)                       # fine window measured directly, matrix form
    qmc = q_needed_matrix(s, Kc + Dp, V' * Lf * V)         # inferred from the coarse lattice, matrix form
    return vcat(sort(eigvals(Symmetric(Kc)); rev = true), sort(eigvals(Symmetric(Kf)); rev = true),
                sort(λc), sort(λf),
                [opnorm(M - I), mis, leak, mis + leak, def,
                 opnorm(M^2 - I), anti, out, tot],
                [dl[end], dl[1], εp, leakage_q_floor(Kf, Lf), s, q_needed(s, kd, ℓ2), q_needed(s, kc, ℓ2)],
                [one_sided_bound(kd, ℓ2, q) for q in qgrid], [one_sided_bound(kc, ℓ2, q) for q in qgrid],
                [θ, km, q_needed(s, km, ℓ2), qmd, qmc], rel)
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

fmt(x, e) = !isfinite(x) ? "   —   " : isfinite(e) ? @sprintf("%.4f(%d)", x, round(Int, e * 1e4)) : @sprintf("%.4f(—)", x)

function main(args)
    cpath, fpath = args[1], args[2]
    opt = Dict("window" => "4,8", "bin" => "10", "drop" => "0", "mref" => "", "qgrid" => "0.1,0.2,0.3", "lag" => "1")
    i = 3
    while i <= length(args)
        key = lstrip(args[i], '-')
        haskey(opt, key) || error("unknown option $(args[i])")
        opt[key] = args[i+1]
        i += 2
    end
    wlev = parse.(Int, split(opt["window"], ","))
    b, drop = parse(Int, opt["bin"]), parse(Int, opt["drop"])
    lag = parse(Int, opt["lag"])
    sref = isempty(opt["mref"]) ? nothing : exp(-lag * parse(Float64, opt["mref"]))
    qgrid = parse.(Float64, split(opt["qgrid"], ","))
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
    ng = length(qgrid)
    per = 4n + 23 + 2ng
    paired = variants == ["fixed", "zero"]
    pix = [4n+2, 4n+4, 4n+5, 4n+9, 4n+12]  # mismatch, bound, defect, total departure, ε₊
    function f(ic, jf)
        q = [window_quantities(meanC(Pc, Mc, ic), meanC(Pf, Mf, jf), sref, qgrid, lag) for (Pf, Mf) in fine]
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
    println("window: stout levels $wlev (ρ = $(kc["rho"])), bins of $b; centered excitations; coarse lag $lag, fine lag $(2lag)")
    for (ia, α) in enumerate(alphas), (iv, var) in enumerate(variants)
        p = ((ia - 1) * length(variants) + (iv - 1)) * per
        println("\nα = $α, $var")
        for (label, m) in (("eig K_c (coarse lag $lag)", n), ("eig K_f (fine lag $(2lag))", n), ("eig L_c", n), ("eig L_f", n))
            println(@sprintf("  %-24s", label), join([fmt(est[p+j], err[p+j]) for j in 1:m], "  "))
            p += m
        end
        for label in ("‖M − I‖ (normalization)", "‖K_f V − V K_c‖", "max leakage^{1/2}",
                      "bound (Proposition 6)", "defect (window_defect)",
                      "|Krein square| ‖M² − I‖", "graph: to anti-aligned", "graph: out of windows",
                      "graph: total departure")
            p += 1
            println(@sprintf("  %-24s %s", label, fmt(est[p], err[p])))
        end
        line(label, vals) = println(@sprintf("  %-36s", label), join(vals, "  "))
        line("one-sided: eig(V*K_fV − K_c)", [fmt(est[p+1], err[p+1]), fmt(est[p+2], err[p+2])])
        line("  ε₊ = max(0, λ_max)", [fmt(est[p+3], err[p+3])])
        line("  q floor from positivity", [fmt(est[p+4], err[p+4])])
        line("  reference s", [fmt(est[p+5], err[p+5])])
        line("  q needed, k = fine window top", [fmt(est[p+6], err[p+6])])
        line("  q needed, k = coarse top + ε₊", [fmt(est[p+7], err[p+7])])
        line("  μ(q), fine k, q = $(join(qgrid, ", "))", [fmt(est[p+7+j], err[p+7+j]) for j in 1:ng])
        line("  μ(q), coarse k", [fmt(est[p+7+ng+j], err[p+7+ng+j]) for j in 1:ng])
        o = p + 7 + 2ng
        line("  angle of leading directions (°)", [@sprintf("%.1f(%.1f)", est[o+1], err[o+1])])
        line("  k = λ_max(K_c + D₊)", [fmt(est[o+2], err[o+2])])
        line("  q needed, k = λ_max(K_c + D₊)", [fmt(est[o+3], err[o+3])])
        line("  q needed, matrix form, fine", [fmt(est[o+4], err[o+4])])
        line("  q needed, matrix form, coarse", [fmt(est[o+5], err[o+5])])
        line("  defect / (1 − s), ε₊ / (1 − s)", [fmt(est[o+6], err[o+6]), fmt(est[o+7], err[o+7])])
    end
    if paired
        println("\nfixed origin − zero momentum (paired errors: same configurations and bins)")
        for (ia, α) in enumerate(alphas)
            p = nvar * per + length(pix) * (ia - 1)
            println("  α = $α   Δ‖K_f V − V K_c‖ ", fmt(est[p+1], err[p+1]), "   Δbound ", fmt(est[p+2], err[p+2]),
                    "   Δdefect ", fmt(est[p+3], err[p+3]), "   Δdeparture ", fmt(est[p+4], err[p+4]),
                    "   Δε₊ ", fmt(est[p+5], err[p+5]))
        end
    end
end

main(ARGS)
