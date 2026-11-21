# Compare a coarse lattice with a blocked fine lattice through an operator window
# (Proposition 6 of docs/notes/spectral_transfer.tex).
#
#   julia --project=lattice lattice/examples/window_transfer.jl coarse.txt fine.txt [options]
#
# coarse.txt comes from examples/glueball.jl, fine.txt from examples/blocked_glueball.jl,
# with the same smearing levels and ρ. Options:
#          --window 2,3      operator levels (1-based) spanning the window
#          --bin 10          configurations per bin, on both lattices
# For each blocking parameter α the script reports the coarse compression K_c and
# leakage L_c at coarse lag 1, the fine ones K_f and L_f at fine lag 2, the alignment V
# induced by the correspondence O ↦ O∘B (the polar part of G_f^{1/2} G_c^{-1/2}), the
# normalization mismatch ‖M − I‖ of its positive part, the in-window mismatch
# ‖K_f V − V K_c‖, and the bound of Proposition 6. Errors: delete-one jackknife over the
# bins of each lattice, the two (independent) variances added.

using LatticeYM, Printf, LinearAlgebra

function read_file(path)
    header = ""
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        isempty(strip(line)) && continue
        startswith(line, "#") ? (header = line) : push!(rows, parse.(Float64, split(line)))
    end
    kv = Dict(split(p, "=")[1] => split(p, "=")[2] for p in split(header[3:end]) if occursin("=", p))
    return kv, permutedims(reduce(hcat, rows))
end

binned_products(O, b) = ((P, M) = correlator_products(O); (bin_means(P, b), bin_means(M, b)))
meanC(Pb, Mb, ix) = correlator_from_means(dropdims(sum(Pb[ix, :, :, :]; dims = 1); dims = 1) ./ length(ix),
                                          vec(sum(Mb[ix, :]; dims = 1)) ./ length(ix))

function quantities(Cc, Cf, w)
    Kc, Lc, Gc = window_compression(Cc[w, w, :], 1)
    Kf, Lf, Gf = window_compression(Cf[w, w, :], 2)
    V, M = induced_alignment(Gf, Gc)
    mis = opnorm(Kf * V - V * Kc)
    λc, λf = eigvals(Symmetric(Lc)), eigvals(Symmetric(Lf))
    leak = sqrt(max(maximum(λc), maximum(λf), 0.0))
    return vcat(sort(eigvals(Symmetric(Kc)); rev = true), sort(eigvals(Symmetric(Kf)); rev = true),
                sort(λc), sort(λf), [opnorm(M - I), mis, leak, mis + leak])
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
    opt = Dict("window" => "2,3", "bin" => "10")
    i = 3
    while i <= length(args)
        opt[lstrip(args[i], '-')] = args[i+1]
        i += 2
    end
    w = parse.(Int, split(opt["window"], ","))
    b = parse(Int, opt["bin"])
    kc, Xc = read_file(cpath)
    kf, Xf = read_file(fpath)
    kc["levels"] == kf["levels"] && kc["rho"] == kf["rho"] ||
        error("the two files must use the same smearing levels and ρ")
    levels = parse.(Int, split(kc["levels"], ","))
    nlev, Tc, Tf = length(levels), parse(Int, kc["T"]), parse(Int, kf["T"])
    Tf == 2Tc || @warn "fine time extent is not twice the coarse one" Tc Tf
    alphas = parse.(Float64, split(kf["alphas"], ","))

    Oc = zeros(size(Xc, 1), nlev, Tc)
    for k in 1:nlev, t in 1:Tc
        Oc[:, k, t] = Xc[:, 1 + (k - 1) * Tc + t]
    end
    Pc, Mc = binned_products(Oc, b)
    n = length(w)
    println("coarse: β = $(kc["beta"]), $(kc["L"])³ × $Tc, $(size(Xc, 1)) configurations; ",
            "fine: β = $(kf["beta"]), $(kf["L"])³ × $Tf blocked to $(kf["blocked"])³, $(size(Xf, 1)) configurations")
    println("window: smearing levels ", levels[w], " (ρ = $(kc["rho"])), bins of $b")
    for (ia, α) in enumerate(alphas)
        Of = zeros(size(Xf, 1), nlev, Tf)
        off = 1 + (ia - 1) * nlev * Tf
        for k in 1:nlev, t in 1:Tf
            Of[:, k, t] = Xf[:, off + (k - 1) * Tf + t]
        end
        Pf, Mf = binned_products(Of, b)
        est, err = two_sample_jackknife((ic, jf) -> quantities(meanC(Pc, Mc, ic), meanC(Pf, Mf, jf), w),
                                        size(Pc, 1), size(Pf, 1))
        p = 0
        println("\nα = $α")
        for (label, m) in (("eig K_c (coarse lag 1)", n), ("eig K_f (fine lag 2)", n),
                           ("eig L_c", n), ("eig L_f", n))
            println(@sprintf("  %-24s", label), join([fmt(est[p+j], err[p+j]) for j in 1:m], "  "))
            p += m
        end
        for label in ("‖M − I‖ (normalization)", "‖K_f V − V K_c‖", "max leakage^{1/2}", "bound (Proposition 6)")
            p += 1
            println(@sprintf("  %-24s %s", label, fmt(est[p], err[p])))
        end
    end
end

main(ARGS)
