# Noise floor of the in-window mismatch of Proposition 6, from one ensemble. Pairs of
# independent block-bootstrap resamples are drawn from the same data, so the underlying
# mismatch is zero by construction; the apparent mismatch ‖K₁V − VK₂‖, with V the induced
# alignment between the two resamples, shows what statistical noise alone produces. The two
# resamples may have different sizes (for example a full coarse ensemble against a pilot).
#
#   julia --project=lattice lattice/examples/window_noise_floor.jl data.txt [options]
#
# Options: --window 4,8   stout levels (values) spanning the window
#          --bin 10       configurations per bootstrap block
#          --size 3000    configurations in each pseudo-ensemble
#          --size2 n      configurations in the second one (default: the same as --size)
#          --pairs 2000   number of resampled pairs
#          --seed 1

using LatticeYM, Printf, LinearAlgebra, Random

function main(args)
    path = args[1]
    opt = Dict("window" => "4,8", "bin" => "10", "size" => "3000", "size2" => "", "pairs" => "2000", "seed" => "1")
    i = 2
    while i <= length(args)
        key = lstrip(args[i], '-')
        haskey(opt, key) || error("unknown option $(args[i])")
        opt[key] = args[i+1]
        i += 2
    end
    _, kv, X, dropped = read_measurements(path)
    dropped && @warn "$path: ignoring an unfinished last line (file still being written?)"
    levels, T = parse.(Int, split(kv["levels"], ",")), parse(Int, kv["T"])
    wlev = parse.(Int, split(opt["window"], ","))
    w = [findfirst(==(l), levels) for l in wlev]
    b, nsize, npairs = parse.(Int, (opt["bin"], opt["size"], opt["pairs"]))
    nsize2 = isempty(opt["size2"]) ? nsize : parse(Int, opt["size2"])
    O = zeros(size(X, 1), length(w), T)
    for (j, k) in enumerate(w), t in 1:T
        O[:, j, t] = X[:, 1 + (k - 1) * T + t]
    end
    P, M = correlator_products(O)
    Pb, Mb = bin_means(P, b), bin_means(M, b)
    nb = size(Pb, 1)
    rng = Xoshiro(parse(Int, opt["seed"]))
    draw(n) = (ndraw = n ÷ b; ix = rand(rng, 1:nb, ndraw);
              correlator_from_means(dropdims(sum(Pb[ix, :, :, :]; dims = 1); dims = 1) ./ ndraw,
                                    vec(sum(Mb[ix, :]; dims = 1)) ./ ndraw))
    mis = zeros(npairs)
    for p in 1:npairs
        K1, _, G1 = window_compression(draw(nsize), 1)
        K2, _, G2 = window_compression(draw(nsize2), 1)
        V, _ = induced_alignment(G1, G2)
        mis[p] = opnorm(K1 * V - V * K2)
    end
    sort!(mis)
    q(x) = mis[clamp(round(Int, x * npairs), 1, npairs)]
    @printf("window %s, %d pairs of resamples of %d and %d configurations (blocks of %d) from %d configurations\n",
            wlev, npairs, nsize, nsize2, b, size(X, 1))
    @printf("apparent mismatch with zero true mismatch: median %.4f, 90%% %.4f, 95%% %.4f, 99%% %.4f\n",
            q(0.5), q(0.9), q(0.95), q(0.99))
end

main(ARGS)
