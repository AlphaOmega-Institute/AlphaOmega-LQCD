using Test, Random, LinearAlgebra
using LatticeYM
using LatticeYM: SUBGROUPS, lmul_su2, su2_near_identity, sample_a0, ONE3

const SYM_TWIST = ((1, 2) => 1, (3, 4) => 1)
const HB_OR = (Heatbath(), Overrelaxation(), Overrelaxation(), Overrelaxation())

# Run `algs` in turn per sweep; return the time series of `f(g)`.
function run_chain(lat, β, algs, ntherm, nmeas, seed; f = mean_plaquette)
    rng = Xoshiro(seed)
    g = hot_start(lat, rng)
    for _ in 1:ntherm, alg in algs
        sweep!(g, β, alg, rng)
    end
    series = map(1:nmeas) do _
        foreach(alg -> sweep!(g, β, alg, rng), algs)
        f(g)
    end
    return series, g
end

# |a − b| within `nσ` combined standard errors.
agree(a, ea, b, eb = 0.0; nσ = 4.5) = abs(a - b) <= nσ * sqrt(ea^2 + eb^2)

@testset "LatticeYM" begin

@testset "SU(3) algebra" begin
    rng = Xoshiro(1)
    for _ in 1:200
        U, V = haar_su3(rng), haar_su3(rng)
        @test su3_deviation(U) < 1e-13
        @test Matrix(U * V) ≈ Matrix(U) * Matrix(V)
        @test Matrix(U') ≈ Matrix(U)'
        @test retr(U, V) ≈ real(tr(Matrix(U) * Matrix(V)))
        @test det(U) ≈ det(Matrix(U))
        a, b = su2_near_identity(0.4, rng)
        i, j = SUBGROUPS[rand(rng, 1:3)]
        @test su3_deviation(lmul_su2(a, b, i, j, U)) < 1e-13
    end
    U = haar_su3(rng)
    @test reunitarize(U) ≈ U
    @test su3_deviation(reunitarize(U + 1e-3 * Mat3(randn(rng, ComplexF64, 3, 3)))) < 1e-13
end

@testset "Haar measure" begin
    rng = Xoshiro(2)
    N = 40_000
    t = [tr(haar_su3(rng)) for _ in 1:N]
    @test abs(sum(t) / N) < 4 / sqrt(N)          # ⟨tr U⟩ = 0 with unit variance
    @test abs(sum(abs2, t) / N - 1) < 0.05       # ⟨|tr U|²⟩ = 1
end

@testset "a₀ sampler, α = $α" for α in (0.0, 0.7, 2.5, 12.0)
    # Compare ⟨a₀⟩ with the exact value from quadrature of √(1−a²) e^{αa}.
    rng = Xoshiro(3)
    xs = [sample_a0(α, rng) for _ in 1:100_000]
    grid = range(-1, 1; length = 20_001)
    w = @. sqrt(max(0, 1 - grid^2)) * exp(α * (grid - 1))
    exact = sum(w .* grid) / sum(w)
    m, e, _ = mean_err(xs)
    @test agree(m, e, exact)
end

@testset "single-plaquette integral" begin
    @test single_plaquette_average(1e-4) ≈ 1e-4 / 18 rtol = 1e-3   # strong coupling: β/18
    @test single_plaquette_average(3.0; n = 128) ≈ single_plaquette_average(3.0; n = 256) atol = 1e-13
    @test 0.99 < single_plaquette_average(400.0) < 1
end

@testset "staple matches the action ($name)" for (name, tw) in (("periodic", ()), ("twisted", SYM_TWIST))
    # This is the check that fails for an inconsistent twist: the change in
    # the total action must equal the local change predicted by the staple.
    rng = Xoshiro(4)
    lat = Lattice((4, 4, 4, 4); twist = tw)
    g = hot_start(lat, rng)
    β = 5.7
    corner = nsites(lat)                                # all coordinates = L
    special = vcat(corner, [lat.fwd[ν, corner] for ν in 1:4], [lat.bwd[ν, corner] for ν in 1:4])
    for x in vcat(special, rand(rng, 1:nsites(lat), 8)), μ in 1:4
        A = staple(g, μ, x)
        Uold, Unew = g.U[μ, x], haar_su3(rng)
        S0 = action(g, β)
        g.U[μ, x] = Unew
        S1 = action(g, β)
        g.U[μ, x] = Uold
        @test S1 - S0 ≈ -β / 3 * (retr(Unew, A) - retr(Uold, A)) atol = 1e-8
    end
end

@testset "gauge invariance ($name)" for (name, tw) in (("periodic", ()), ("twisted", SYM_TWIST))
    rng = Xoshiro(5)
    lat = Lattice((4, 4, 4, 6); twist = tw)
    g = hot_start(lat, rng)
    for _ in 1:5
        sweep!(g, 5.7, Heatbath(), rng)                 # move away from random links
    end
    h = gauge_transform!(copy(g), [haar_su3(rng) for _ in 1:nsites(lat)])
    @test action(h, 5.7) ≈ action(g, 5.7) rtol = 1e-12
    @test plane_plaquettes(h) ≈ plane_plaquettes(g) atol = 1e-12
    @test polyakov_loop(h) ≈ polyakov_loop(g) atol = 1e-12
    for (μ, ν) in ((1, 2), (3, 4), (2, 4))
        @test wilson_loop(h, μ, ν, 2, 3) ≈ wilson_loop(g, μ, ν, 2, 3) atol = 1e-12
        @test wilson_loop(g, μ, ν, 1, 1) ≈ plane_plaquettes(g)[μ, ν] atol = 1e-12
    end
end

@testset "overrelaxation conserves the action" begin
    rng = Xoshiro(6)
    for tw in ((), SYM_TWIST)
        lat = Lattice((4, 4, 4, 4); twist = tw)
        g = hot_start(lat, rng)
        for _ in 1:5
            sweep!(g, 5.7, Heatbath(), rng)
        end
        S0 = action(g, 5.7)
        sweep!(g, 5.7, Overrelaxation(), rng)
        @test action(g, 5.7) ≈ S0 rtol = 1e-11
    end
end

@testset "reproducible with a fixed seed" begin
    # 8⁴ is large enough to take the multithreaded path when julia runs with -t > 1.
    lat = Lattice((8, 8, 8, 8))
    s1, _ = run_chain(lat, 5.7, HB_OR, 0, 2, 7)
    s2, _ = run_chain(lat, 5.7, HB_OR, 0, 2, 7)
    @test s1 == s2
end

# The exact 2D plaquette. Finite-volume corrections on a 16×16 torus are of
# order u^256 and invisible.
@testset "2D exact plaquette, β = $β, $name" for β in (1.0, 4.0, 8.0),
        (name, algs) in (("heatbath+OR", HB_OR), ("Metropolis", (Metropolis(ε = 0.5, nhit = 8),)))
    lat = Lattice((16, 16))
    series, _ = run_chain(lat, β, algs, 300, 4000, 10 + round(Int, β))
    m, e, _ = mean_err(series)
    @test agree(m, e, single_plaquette_average(β))
end

@testset "2D twisted: same exact plaquette, corner = bulk" begin
    lat = Lattice((16, 16); twist = ((1, 2) => 1,))
    β = 4.0
    series, _ = run_chain(lat, β, HB_OR, 300, 4000, 20; f = g -> corner_and_bulk(g, 1, 2))
    mc, ec, _ = mean_err(first.(series))
    mb, eb, _ = mean_err(last.(series))
    exact = single_plaquette_average(β)
    @test agree(mb, eb, exact)
    @test agree(mc, ec, exact)
end

@testset "4D strong coupling, β = 1" begin
    # ⟨P⟩ = u(β) + O(u⁵) with u ≈ 0.06, so the correction is ~1e-6.
    lat = Lattice((4, 4, 4, 4))
    series, _ = run_chain(lat, 1.0, HB_OR, 100, 2000, 30)
    m, e, _ = mean_err(series)
    @test agree(m, e, single_plaquette_average(1.0))
end

@testset "4D heatbath and Metropolis agree, β = 5.7" begin
    lat = Lattice((4, 4, 4, 4))
    s1, _ = run_chain(lat, 5.7, HB_OR, 200, 1500, 40)
    s2, _ = run_chain(lat, 5.7, (Metropolis(ε = 0.35, nhit = 6),), 400, 3000, 41)
    m1, e1, _ = mean_err(s1)
    m2, e2, _ = mean_err(s2)
    @test agree(m1, e1, m2, e2)
end

@testset "4D twist: corner plaquettes equal the bulk" begin
    # A twist that gives both orientations of a plaquette the same phase fails this badly.
    lat = Lattice((4, 4, 4, 4); twist = SYM_TWIST)
    series, _ = run_chain(lat, 5.7, HB_OR, 200, 2000, 50;
                          f = g -> (corner_and_bulk(g, 1, 2)..., corner_and_bulk(g, 3, 4)...))
    for k in (1, 3)
        mc, ec, _ = mean_err(getindex.(series, k))
        mb, eb, _ = mean_err(getindex.(series, k + 1))
        @test agree(mc, ec, mb, eb)
    end
end

# ---------------------------------------------------------------------------
# Smearing, glueball operators, Wilson loops, correlator analysis

# a thermalized-looking field on a small 4D lattice
function rough_field(dims, seed; β = 5.7, n = 6)
    rng = Xoshiro(seed)
    g = hot_start(Lattice(dims), rng)
    for _ in 1:n
        sweep!(g, β, Heatbath(), rng)
    end
    return g, rng
end

@testset "expm3 matches the dense exponential" begin
    rng = Xoshiro(60)
    for s in (0.01, 0.3, 1.0, 4.0, 12.0)
        for _ in 1:20
            M = randn(rng, ComplexF64, 3, 3)
            H = (M + M') / 2
            H -= tr(H) / 3 * I
            X = Mat3(1im * s * H / norm(H))
            E = LatticeYM.expm3(X)
            @test Matrix(E) ≈ exp(Matrix(X)) atol = 1e-12
            @test su3_deviation(E) < 1e-12
        end
    end
end

@testset "stout smearing: SU(3), time links untouched, unit and pure-gauge fields fixed" begin
    g, rng = rough_field((4, 4, 4, 6), 61)
    h = stout_spatial(g, 0.1, 5)
    @test maximum(su3_deviation, h.U[1:3, :]) < 1e-12
    @test h.U[4, :] == g.U[4, :]
    @test all(i -> isapprox(stout_spatial(g, 0.0).U[i], g.U[i]; atol = 1e-14), eachindex(g.U))
    c = cold_start(Lattice((4, 4, 4, 6)))
    @test maximum(su3_deviation, stout_spatial(c, 0.2, 3).U) < 1e-14
    @test all(x -> isapprox(x, ONE3; atol = 1e-14), stout_spatial(c, 0.2, 3).U)
    pure = gauge_transform!(copy(c), [haar_su3(rng) for _ in 1:nsites(c.lat)])
    @test all(≈(1.0; atol = 1e-12), slice_plaquettes(stout_spatial(pure, 0.15, 4)))
end

@testset "stout smearing commutes with gauge transformations" begin
    g, rng = rough_field((4, 4, 4, 6), 62)
    G = [haar_su3(rng) for _ in 1:nsites(g.lat)]
    a = gauge_transform!(stout_spatial(g, 0.12, 3), G)
    b = stout_spatial(gauge_transform!(copy(g), G), 0.12, 3)
    @test maximum(maximum(abs, (a.U[i] - b.U[i]).e) for i in eachindex(a.U)) < 1e-12
end

@testset "stout smearing smooths: spatial plaquettes increase" begin
    g, _ = rough_field((4, 4, 4, 6), 63)
    p = [sum(slice_plaquettes(stout_spatial(g, 0.1, n))) for n in (0, 1, 2, 4, 8)]
    @test issorted(p)
end

@testset "glueball operators and Wilson loops" begin
    g, rng = rough_field((4, 4, 4, 6), 64)
    # slice average = mean of the spatial planes
    P = plane_plaquettes(g)
    @test sum(slice_plaquettes(g)) / 6 ≈ (P[1, 2] + P[1, 3] + P[2, 3]) / 3 atol = 1e-12
    levels = [0, 2, 5]
    O = glueball_operators(g, 0.1, levels)
    for (k, n) in enumerate(levels)
        @test O[k, :] ≈ slice_plaquettes(stout_spatial(g, 0.1, n)) atol = 1e-13
    end
    G = [haar_su3(rng) for _ in 1:nsites(g.lat)]
    gt = gauge_transform!(copy(g), G)
    @test glueball_operators(gt, 0.1, levels) ≈ O atol = 1e-12
    W = wilson_loops_rt(g, 3, 4)
    for R in 1:3, T in 1:4
        @test W[R, T] ≈ sum(wilson_loop(g, μ, 4, R, T) for μ in 1:3) / 3 atol = 1e-12
    end
    @test W[1, 1] ≈ (P[1, 4] + P[2, 4] + P[3, 4]) / 3 atol = 1e-12
    @test wilson_loops_rt(gt, 3, 4) ≈ W atol = 1e-12
end

@testset "correlator analysis on exact and synthetic data" begin
    rng = Xoshiro(65)
    T = 24
    # effective mass of an exact periodic correlator
    for m in (0.3, 0.77, 1.4)
        c = [cosh(m * (t - T / 2)) for t in 0:T-1]
        @test all(≈(m; atol = 1e-9), effective_masses(c))
    end
    # GEVP recovers both energies of an exact two-state, two-operator correlator
    E = [0.8, 1.6]; Z = [1.0 0.5; 0.7 -0.9]
    C = zeros(2, 2, T)
    for t in 0:T-1
        C[:, :, t+1] = Z * Diagonal(exp.(-E .* t)) * Z'
    end
    for (t0, t) in ((0, 1), (1, 3), (2, 5))
        λ, V = gevp(C, t0, t)
        @test λ ≈ exp.(-E .* (t - t0)) rtol = 1e-10
        @test V' * C[:, :, t0+1] * V ≈ I atol = 1e-10
    end
    # correlator_matrix against a direct implementation
    O = randn(rng, 5, 3, 8)
    Cm = correlator_matrix(O)
    m = vec(sum(O; dims = (1, 3))) ./ (5 * 8)
    for k in 1:3, l in 1:3, t in 0:7
        direct = sum(O[c, k, mod1(t0 + t, 8)] * O[c, l, t0] for c in 1:5, t0 in 1:8) / 40 - m[k] * m[l]
        directT = sum(O[c, l, mod1(t0 + t, 8)] * O[c, k, t0] for c in 1:5, t0 in 1:8) / 40 - m[k] * m[l]
        @test Cm[k, l, t+1] ≈ (direct + directT) / 2 atol = 1e-12
    end
    # products binned before the subtraction reproduce the unbinned correlator
    P, M = correlator_products(O)
    Pm = dropdims(sum(P; dims = 1); dims = 1) ./ 5
    Mm = vec(sum(M; dims = 1)) ./ 5
    @test correlator_from_means(Pm, Mm) ≈ Cm atol = 1e-12
    # transfer residual: exact for a finite positive transfer matrix and a trial vector
    for _ in 1:20
        Q = Matrix(qr(randn(rng, 6, 6)).Q)
        Tm = Q * Diagonal(rand(rng, 6)) * Q'                  # positive, norm < 1
        ψ = randn(rng, 6); ψ /= norm(ψ)
        c = [ψ' * Tm^t * ψ for t in 0:3] .* 2.7              # arbitrary normalization
        exact = norm((I - ψ * ψ') * Tm * ψ)
        @test transfer_residual(c) ≈ exact atol = 1e-12
    end
    # jackknife error of a mean equals the standard error; binning averages blocks
    x = randn(rng, 50)
    est, err = jackknife(v -> sum(v) / length(v), x)
    @test est ≈ sum(x) / 50
    @test err ≈ sqrt(sum(abs2, x .- est) / (49 * 50)) rtol = 1e-12
    @test bin_means(collect(1.0:10.0), 3) ≈ [2.0, 5.0, 8.0]
end

end
