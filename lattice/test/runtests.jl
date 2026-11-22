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


@testset "su3_project: SU(3), fixed points, covariance at both ends" begin
    rng = Xoshiro(70)
    for _ in 1:50
        A = Mat3(randn(rng, ComplexF64, 3, 3))
        P = su3_project(A)
        @test su3_deviation(P) < 1e-12
        G, H = haar_su3(rng), haar_su3(rng)
        @test maximum(abs, (su3_project(G * A * H') - G * P * H').e) < 1e-11
        U = haar_su3(rng)
        @test maximum(abs, (su3_project(U) - U).e) < 1e-12
        @test maximum(abs, (su3_project(2.5 * U) - U).e) < 1e-12
    end
    # row Gram–Schmidt is not covariant under left multiplication, which is why blocking does not use it
    A = Mat3(randn(rng, ComplexF64, 3, 3)); G = haar_su3(rng)
    @test maximum(abs, (reunitarize(G * A) - G * reunitarize(A)).e) > 1e-3
end

@testset "block_spatial: geometry, α = 0, unit field, gauge covariance" begin
    g, rng = rough_field((8, 8, 8, 6), 71)
    b = block_spatial(g, 0.3)
    @test b.lat.dims == (4, 4, 4, 6)
    @test maximum(su3_deviation, b.U) < 1e-12
    b0 = block_spatial(g, 0.0)
    li = LinearIndices(g.lat.dims)
    for X in (1, 17, 40, nsites(b0.lat))
        c = b0.lat.coords[:, X]
        x = li[2c[1] - 1, 2c[2] - 1, 2c[3] - 1, c[4]]
        for μ in 1:3
            @test maximum(abs, (b0.U[μ, X] - g.U[μ, x] * g.U[μ, g.lat.fwd[μ, x]]).e) < 1e-12
        end
        @test b0.U[4, X] == g.U[4, x]
    end
    c = cold_start(Lattice((8, 8, 8, 6)))
    @test all(u -> maximum(abs, (u - ONE3).e) < 1e-12, block_spatial(c, 0.4).U)
    G = [haar_su3(rng) for _ in 1:nsites(g.lat)]
    gb = block_spatial(gauge_transform!(copy(g), G), 0.3)
    Gc = [G[li[2b.lat.coords[1, X] - 1, 2b.lat.coords[2, X] - 1, 2b.lat.coords[3, X] - 1, b.lat.coords[4, X]]]
          for X in 1:nsites(b.lat)]
    bg = gauge_transform!(copy(b), Gc)
    @test maximum(maximum(abs, (gb.U[i] - bg.U[i]).e) for i in eachindex(gb.U)) < 1e-11
    @test glueball_operators(gb, 0.1, [0, 3]) ≈ glueball_operators(b, 0.1, [0, 3]) atol = 1e-11
end

@testset "operator windows: compressions, leakage, alignment, Proposition 6" begin
    rng = Xoshiro(72)
    psd(n) = (A = randn(rng, n, n); M = A * A'; Symmetric(M / opnorm(M)))
    isq(G) = (F = eigen(Symmetric(G)); F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors')
    for _ in 1:30
        nf, nc, k, T = 20, 9, rand(rng, 1:3), 12
        Sf, Sc = psd(nf), psd(nc)
        F, Cs = randn(rng, nf, k), randn(rng, nc, k)
        Cf = cat((F' * Sf^t * F for t in 0:T-1)...; dims = 3)
        Cc = cat((Cs' * Sc^t * Cs for t in 0:T-1)...; dims = 3)
        Kf, Lf, Gf = window_compression(Cf, 1)
        Kc, Lc, Gc = window_compression(Cc, 1)
        Uf, Uc = F * isq(F' * F), Cs * isq(Cs' * Cs)
        @test Kf ≈ Uf' * Sf * Uf atol = 1e-10
        @test Lf ≈ Uf' * Sf * (I - Uf * Uf') * Sf * Uf atol = 1e-10
        V, M = induced_alignment(Gf, Gc)
        @test V' * V ≈ I atol = 1e-10
        @test V * M ≈ sqrt(Symmetric(Gf)) * isq(Gc) atol = 1e-10
        W = Uf * V * Uc'
        @test opnorm(Sf * W - W * Sc) <= window_defect_bound(Kf, Lf, Kc, Lc, V) + 1e-10
        # ℓ = 2 uses lags 0, 2, 4
        K2, L2, _ = window_compression(Cf, 2)
        @test K2 ≈ Uf' * Sf^2 * Uf atol = 1e-10
    end
    _, M = induced_alignment([2.0 0.3; 0.3 1.0], [2.0 0.3; 0.3 1.0])
    @test M ≈ I atol = 1e-12
end

@testset "ground_state_fractions: exact at every t0" begin
    E = [0.8, 1.7]; Z = [1.0 0.6; 0.4 -0.9]                  # Z[k, n]: operator k, state n
    for periodic in (false, true), t0 in (0, 1, 2)
        T = 24
        C = zeros(2, 2, T)
        for t in 0:T-1
            w = periodic ? exp.(-E .* t) .+ exp.(-E .* (T - t)) : exp.(-E .* t)
            C[:, :, t+1] = Z * Diagonal(w) * Z'
        end
        f, E0 = ground_state_fractions(C, [1, 2]; t0, periodic)
        @test E0 ≈ E[1] rtol = 1e-8
        w0 = periodic ? 1 + exp(-E[1] * T) : 1.0
        @test f ≈ [Z[k, 1]^2 * w0 / C[k, k, 1] for k in 1:2] rtol = 1e-8
    end
    # a single exponential with overlap one gives f = 1 at any t0
    C1 = reshape([3.0 * exp(-0.9t) for t in 0:15], 1, 1, 16)
    for t0 in 0:3
        @test ground_state_fractions(C1, [1]; t0, periodic = false)[1][1] ≈ 1 rtol = 1e-10
    end
end


@testset "su3_project fallback for singular path sums" begin
    rng = Xoshiro(73)
    M = randn(rng, ComplexF64, 3, 2) * randn(rng, ComplexF64, 2, 3)       # rank two
    U = haar_su3(rng)
    @test su3_project(Mat3(M), U) == U
    @test_throws ArgumentError su3_project(Mat3(M))
    A = Mat3(randn(rng, ComplexF64, 3, 3))
    @test su3_project(A, U) ≈ su3_project(A)                               # fallback unused
end

@testset "su3_project at the determinant branch cut" begin
    rng = Xoshiro(76)
    on(X) = opnorm(Matrix(X))
    V, F = haar_su3(rng), haar_su3(rng)
    B = randn(rng, ComplexF64, 3, 3)
    P = B' * B + I                                       # positive definite
    A = Mat3(-Matrix(V) * P)                             # polar part −V: det = −1, on the cut
    n0 = LatticeYM.SU3_FALLBACKS[]
    @test su3_project(A, F) == F
    @test LatticeYM.SU3_FALLBACKS[] == n0 + 1
    @test_throws ArgumentError su3_project(A)
    # with no band (cut_tol < 0, the principal root always) roundoff picks either cube root,
    # so gauge-equivalent inputs project to links differing by a center phase, norm √3
    d = map(1:200) do _
        G, H = haar_su3(rng), haar_su3(rng)
        on(G' * su3_project(G * A * H'; cut_tol = -1) * H - su3_project(A; cut_tol = -1))
    end
    @test any(x -> abs(x - sqrt(3)) < 1e-9, d) && all(x -> x < 1e-9 || abs(x - sqrt(3)) < 1e-9, d)
    # with the band every gauge-equivalent input gets the transformed fallback
    for _ in 1:200
        G, H = haar_su3(rng), haar_su3(rng)
        @test on(su3_project(G * A * H', G * F * H') - G * F * H') < 1e-12
    end
    # just outside the band the projection is used, and is covariant
    A2 = Mat3(cis((π - 1e-6) / 3) * Matrix(V) * P)
    n1 = LatticeYM.SU3_FALLBACKS[]
    @test on(su3_project(A2, F) - V) < 1e-12
    G, H = haar_su3(rng), haar_su3(rng)
    @test on(su3_project(G * A2 * H', G * F * H') - G * V * H') < 1e-11
    @test LatticeYM.SU3_FALLBACKS[] == n1

    # a lattice configuration with a blocked path sum on the cut, at α = 0.6: unit links except
    # U_1(x0) = diag(e^{iφ}, e^{0.1i}, e^{−i(φ+0.1)}) and U_2(x0) = U_3(x0) = e^{2πi/3}, so the path
    # sum of the coarse link (1, x0) is U_1(x0) + 2α e^{iπ/3}; φ solves det(polar part) = −1
    α = 0.6
    w = 2α * cis(π / 3)
    dp(φ) = prod(z -> (z + w) / abs(z + w), (cis(φ), cis(0.1), cis(-φ - 0.1)))
    a, b = 1.9, 2.1
    @assert imag(dp(a)) * imag(dp(b)) < 0 && real(dp((a + b) / 2)) < 0
    for _ in 1:100
        m = (a + b) / 2
        imag(dp(a)) * imag(dp(m)) <= 0 ? (b = m) : (a = m)
    end
    φ = (a + b) / 2
    @test abs(dp(φ) + 1) < 1e-14
    g = cold_start(Lattice((8, 8, 8, 4)))
    li = LinearIndices(g.lat.dims)
    x0 = li[1, 1, 1, 1]
    g.U[1, x0] = Mat3(Matrix(Diagonal([cis(φ), cis(0.1), cis(-φ - 0.1)])))
    g.U[2, x0] = g.U[3, x0] = cis(2π / 3) * ONE3
    n2 = LatticeYM.SU3_FALLBACKS[]
    bl = block_spatial(g, α)
    @test LatticeYM.SU3_FALLBACKS[] == n2 + 1                 # this one link, and no other
    @test bl.lat.coords[:, 1] == [1, 1, 1, 1] && bl.U[1, 1] == g.U[1, x0]   # its straight path
    for _ in 1:10
        G = [haar_su3(rng) for _ in 1:nsites(g.lat)]
        gb = block_spatial(gauge_transform!(copy(g), G), α)
        Gc = [G[li[2bl.lat.coords[1, X] - 1, 2bl.lat.coords[2, X] - 1, 2bl.lat.coords[3, X] - 1, bl.lat.coords[4, X]]]
              for X in 1:nsites(bl.lat)]
        bg = gauge_transform!(copy(bl), Gc)
        @test maximum(maximum(abs, (gb.U[i] - bg.U[i]).e) for i in eachindex(gb.U)) < 1e-11
    end
end

@testset "measurement files: complete lines only, every row checked, prefix hash" begin
    path = joinpath(mktempdir(), "m.txt")
    head = "# beta=6.0 L=4 levels=1,2"
    body = head * "\n1.0 2.0 3.0\n4.0 5.0 6.0\n"
    write(path, body)
    h, kv, X, dropped = read_measurements(path)
    @test h == head && kv["beta"] == "6.0" && kv["levels"] == "1,2"
    @test X == [1.0 2.0 3.0; 4.0 5.0 6.0] && !dropped
    # a writer caught mid-line: an exponent cut short, or a line that happens to have the
    # full width or to end in a valid number, is not read
    for tail in ("7.0 8.0 1.2e-", "7.0 8.0 1.2", "7.0 8.0 9.0 1.2", "7")
        write(path, body * tail)
        _, _, X, dropped = read_measurements(path)
        @test size(X) == (2, 3) && dropped
    end
    # every complete row is checked, and two runs joined together are refused
    for bad in ("4.0 5.0e- 6.0\n7.0 8.0 9.0\n", "4.0 5.0\n7.0 8.0 9.0\n", head * "\n7.0 8.0 9.0\n", "4.0 x 6.0\n")
        write(path, head * "\n1.0 2.0 3.0\n" * bad)
        @test_throws ErrorException read_measurements(path)
    end
    write(path, "1.0 2.0\n")
    @test_throws ErrorException read_measurements(path)          # no header
    # the FNV-1a prefix hash: published test vectors, pieces, and file prefixes
    @test fnv1a("") == 0xcbf29ce484222325 && fnv1a("a") == 0xaf63dc4c8601ec8c && fnv1a("foobar") == 0x85944171f73967e8
    write(path, body)
    bytes = read(path)
    @test file_prefix_hash(path, length(bytes)) == fnv1a(bytes) == fnv1a(body)
    @test file_prefix_hash(path, 10) == fnv1a(bytes[1:10])
    @test fnv1a(bytes[11:end], fnv1a(bytes[1:10])) == fnv1a(bytes)
    @test_throws ErrorException file_prefix_hash(path, length(bytes) + 1)
end

@testset "blocking origins and the zero-momentum projection" begin
    g, rng = rough_field((8, 8, 8, 6), 74)
    li = LinearIndices(g.lat.dims)
    b = block_spatial(g, 0.0; origin = (1, 0, 1))
    X = 23
    c = b.lat.coords[:, X]
    x = li[2c[1], 2c[2] - 1, 2c[3], c[4]]
    @test maximum(abs, (b.U[2, X] - g.U[2, x] * g.U[2, g.lat.fwd[2, x]]).e) < 1e-12
    G = [haar_su3(rng) for _ in 1:nsites(g.lat)]
    gb = block_spatial(gauge_transform!(copy(g), G), 0.3; origin = (0, 1, 1))
    bb = block_spatial(g, 0.3; origin = (0, 1, 1))
    Gc = [G[li[2bb.lat.coords[1, X] - 1, 2bb.lat.coords[2, X], 2bb.lat.coords[3, X], bb.lat.coords[4, X]]]
          for X in 1:nsites(bb.lat)]
    bg = gauge_transform!(copy(bb), Gc)
    @test maximum(maximum(abs, (gb.U[i] - bg.U[i]).e) for i in eachindex(gb.U)) < 1e-11
    # translating the fine field by one site changes the fixed-origin operator, but not the
    # average over the eight origins
    shifted = GaugeField(g.lat, [g.U[μ, g.lat.bwd[1, x]] for μ in 1:4, x in 1:nsites(g.lat)])
    origins = vec([(a, b, c) for a in 0:1, b in 0:1, c in 0:1])
    zero_mom(h) = sum(glueball_operators(block_spatial(h, 0.3; origin = o), 0.1, [0, 2]) for o in origins) ./ 8
    fixed(h) = glueball_operators(block_spatial(h, 0.3), 0.1, [0, 2])
    @test zero_mom(shifted) ≈ zero_mom(g) atol = 1e-12
    @test !isapprox(fixed(shifted), fixed(g); atol = 1e-6)
end

@testset "window_defect equals the defect of the window map" begin
    rng = Xoshiro(75)
    psd(n) = (A = randn(rng, n, n); M = A * A'; Symmetric(M / opnorm(M)))
    isq(G) = (F = eigen(Symmetric(G)); F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors')
    for _ in 1:40
        nf, nc, k, T = 18, 10, rand(rng, 1:4), 10
        Sf, Sc = psd(nf), psd(nc)
        F, Cs = randn(rng, nf, k), randn(rng, nc, k)
        Cf = cat((F' * Sf^t * F for t in 0:T-1)...; dims = 3)
        Cc = cat((Cs' * Sc^t * Cs for t in 0:T-1)...; dims = 3)
        Kf, Lf, Gf = window_compression(Cf, 1)
        Kc, Lc, Gc = window_compression(Cc, 1)
        V, _ = induced_alignment(Gf, Gc)
        W = F * isq(F' * F) * V * isq(Cs' * Cs) * Cs'
        exact = opnorm(Sf * W - W * Sc)
        @test window_defect(Kf, Lf, Kc, Lc, V) ≈ exact rtol = 1e-8
        @test exact <= window_defect_bound(Kf, Lf, Kc, Lc, V) + 1e-10
        Uc = Cs * isq(Cs' * Cs)                                    # coarse-window inputs only
        @test window_defect_inputs(Kf, Lf, Kc, V) ≈ opnorm((Sf * W - W * Sc) * Uc) rtol = 1e-7
        @test window_defect_inputs(Kf, Lf, Kc, V) <= exact + 1e-10
    end
end

@testset "one-sided gap bound, positivity floor on q, q needed" begin
    rng = Xoshiro(78)
    for _ in 1:200
        n, k = 12, rand(rng, 1:3)
        Ω = normalize(randn(rng, n))
        Pperp = I - Ω * Ω'
        A = Pperp * randn(rng, n, n) * Pperp
        R = A * A'
        R = R / opnorm(R) * rand(rng, 0.2:0.01:0.95)            # spectrum on Ω⊥ below 1
        T = Ω * Ω' + R                                          # 0 ≤ T ≤ 1, TΩ = Ω
        U = Matrix(qr(Pperp * randn(rng, n, k)).Q)[:, 1:k]      # window ⊥ Ω
        Q = Pperp - U * U'
        K = U' * T * U
        L = U' * T^2 * U - K^2                                  # = (QTU)*(QTU)
        kk, ℓ2, q = eigmax(Symmetric(K)), opnorm(Q * T * U)^2, opnorm(Q * T * Q)
        @test ℓ2 ≈ eigmax(Symmetric((L + L') / 2)) rtol = 1e-8 atol = 1e-12
        λ2 = opnorm(R)
        @test λ2 <= one_sided_bound(kk, ℓ2, q) + 1e-12
        @test leakage_q_floor(K, L) <= q + 1e-10                # L ≤ qK
        s = kk + rand(rng) * (1 - kk)
        qn = q_needed(s, kk, ℓ2)
        qn >= 0 && @test one_sided_bound(kk, ℓ2, qn) ≈ s rtol = 1e-10
    end
    @test isnan(q_needed(0.3, 0.4, 0.01))
    @test q_needed(0.3, 0.3, 0.0) == 0.3 && one_sided_bound(0.3, 0.0, 0.3) ≈ 0.3   # k = s, no leakage
    @test q_needed(0.3, 0.2, -1e-4) == 0.3 && isnan(q_needed(0.3, 0.3, 1e-6))
    @test isnan(leakage_q_floor([1.0 0; 0 -0.1], [0.1 0; 0 0.1]))
end

@testset "operator-valued one-sided step and q_needed_matrix" begin
    rng = Xoshiro(80)
    for _ in 1:200
        n, k = 14, rand(rng, 1:4)
        Ω = normalize(randn(rng, n))
        Pp = I - Ω * Ω'
        A = Pp * randn(rng, n, n) * Pp
        R = A * A'
        R = R / opnorm(R) * 0.9
        T = Ω * Ω' + R
        U = Matrix(qr(Pp * randn(rng, n, k)).Q)[:, 1:k]
        Qb = nullspace([Ω U]')
        T2 = T^2
        C = U' * T2 * U
        L = U' * T2 * (Qb * Qb') * T2 * U
        q = opnorm(Qb' * T2 * Qb)
        X = [U Qb]
        γ = 10.0^(4rand(rng) - 2)                 # T² ≤ (C + L/γ) ⊕ (q + γ) I on the vacuum complement
        Bm = [C + L / γ zeros(k, size(Qb, 2)); zeros(size(Qb, 2), k) (q + γ) * I]
        @test eigmin(Symmetric(Bm - X' * T2 * X)) > -1e-10
        s = q + 0.05 + rand(rng) * 0.5
        qn = q_needed_matrix(s, C, L)
        if isfinite(qn) && q <= qn                 # the test passes, so T² ≤ s
            @test opnorm(X' * T2 * X) <= s + 1e-10
        end
        if isfinite(qn)                            # the matrix form never needs a smaller q than the scalar form
            @test qn >= q_needed(s, eigmax(Symmetric(C)), eigmax(Symmetric(L))) - 1e-10
        end
    end
    @test q_needed_matrix(0.4, fill(0.3, 1, 1), fill(0.01, 1, 1)) ≈ q_needed(0.4, 0.3, 0.01) rtol = 1e-10
    @test isnan(q_needed_matrix(0.3, fill(0.35, 1, 1), fill(0.01, 1, 1)))
    # λ_max(C) = s: passes when the leakage vanishes along the saturating direction, fails otherwise
    @test q_needed_matrix(0.4, [0.4 0; 0 0.2], [0 0; 0 0.01]) ≈ 0.35 rtol = 1e-9
    @test isnan(q_needed_matrix(0.4, [0.4 0; 0 0.2], [0.01 0; 0 0.01]))
end

@testset "correlated cosh fit" begin
    rng = Xoshiro(79)
    for m in (0.4, 0.66, 1.0), T in (16, 32)
        c = [3.0 * cosh(m * (t - T / 2)) for t in 0:T-1]
        ts = 1:T÷4
        X = randn(rng, length(ts), length(ts))
        cov = (X * X' + I) .* 1e-6
        mf, A, χ = cosh_fit(c, cov, ts, T)
        @test mf ≈ m rtol = 1e-7
        @test A ≈ 3.0 rtol = 1e-6
        @test χ < 1e-8
    end
    # with noise drawn from the covariance, χ² per degree of freedom is of order one on average
    T, m, ts = 32, 0.66, 2:8
    c0 = [cosh(m * (t - T / 2)) for t in 0:T-1]
    cov = Diagonal((0.02 .* c0[ts .+ 1]) .^ 2) |> Matrix
    χs, ms = Float64[], Float64[]
    for _ in 1:200
        c = copy(c0)
        c[ts .+ 1] .+= sqrt.(diag(cov)) .* randn(rng, length(ts))
        mf, _, χ = cosh_fit(c, cov, ts, T)
        push!(χs, χ); push!(ms, mf)
    end
    @test 0.7 < sum(χs) / 200 / (length(ts) - 2) < 1.3
    @test abs(sum(ms) / 200 - m) < 0.01
end

@testset "cosh_mass over one step and over longer intervals" begin
    for m in (0.3, 0.7, 1.4), T in (16, 32), lag in (1, 2, 3)
        c = [cosh(m * (t - T / 2)) for t in 0:T-1]
        for t in 0:T÷2-lag
            @test cosh_mass(c[t+1] / c[t+lag+1], t, T; lag) ≈ m rtol = 1e-10
        end
    end
    @test isnan(cosh_mass(0.9, 1, 16; lag = 2))
end

@testset "the window map in the doubled space: reflection, graphs, normalization" begin
    rng = Xoshiro(77)
    psd(n) = (A = randn(rng, n, n); M = A * A'; Symmetric(M / opnorm(M)))
    isq(G) = (F = eigen(Symmetric(G)); F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors')
    for _ in 1:40
        nf, nc, k, T = 18, 10, rand(rng, 1:4), 10
        Sf, Sc = psd(nf), psd(nc)
        F, Cs = randn(rng, nf, k), randn(rng, nc, k)
        Cf = cat((F' * Sf^t * F for t in 0:T-1)...; dims = 3)
        Cc = cat((Cs' * Sc^t * Cs for t in 0:T-1)...; dims = 3)
        Kf, Lf, Gf = window_compression(Cf, 1)
        Kc, Lc, Gc = window_compression(Cc, 1)
        V, M = induced_alignment(Gf, Gc)
        Uf, Uc = F * isq(F' * F), Cs * isq(Cs' * Cs)          # orthonormal window frames
        W = Uf * V * Uc'
        # coarse ⊕ fine: doubled evolution and the polar reflection
        S = [Sc zeros(nc, nf); zeros(nf, nc) Sf]
        R = [zeros(nc, nc) W'; W zeros(nf, nf)]
        @test opnorm(S * R - R * S) ≈ opnorm(Sf * W - W * Sc) rtol = 1e-10
        Gp, Gm = [Uc; Uf * V] / sqrt(2), [Uc; -Uf * V] / sqrt(2)  # aligned and anti-aligned graphs
        @test R * Gp ≈ Gp && R * Gm ≈ -Gm
        Pp = Gp * Gp'
        Pw = [Uc * Uc' zeros(nc, nf); zeros(nf, nc) Uf * Uf']   # both windows = P₊ + P₋
        @test Pw ≈ Pp + Gm * Gm'
        anti, out, tot = graph_departure(Kf, Lf, Kc, Lc, V)
        @test anti ≈ opnorm(Gm' * S * Gp) rtol = 1e-8
        @test anti ≈ opnorm(Kf * V - V * Kc) / 2 rtol = 1e-8
        @test out ≈ opnorm((I - Pw) * S * Gp) rtol = 1e-7 atol = 1e-9
        @test tot ≈ opnorm((I - Pp) * S * Gp) rtol = 1e-8
        # raw correspondence A = V M: largest change of length, and largest absolute Krein square
        # ‖x‖² − ‖Ax‖² of its graph {(x, Ax)} for unit x
        A = V * M
        x = svd(M).V
        @test opnorm(M - I) ≈ maximum(abs(norm(A * x[:, j]) - 1) for j in 1:k) rtol = 1e-10
        @test opnorm(M^2 - I) ≈ opnorm(A' * A - I) rtol = 1e-10
    end
    # complex Hermitian windows: the formalism allows them (lattice correlators are real)
    hpsd(n) = (A = randn(rng, ComplexF64, n, n); M = A * A'; Matrix(Hermitian(M / opnorm(M))))
    frame(n, k) = Matrix(qr(randn(rng, ComplexF64, n, k)).Q)[:, 1:k]
    for _ in 1:20
        nf, nc, k = 7, 5, 2
        Sf, Sc = hpsd(nf), hpsd(nc)
        Uf, Uc, V = frame(nf, k), frame(nc, k), frame(k, k)
        Kf, Kc = Uf' * Sf * Uf, Uc' * Sc * Uc
        Lf, Lc = Uf' * Sf^2 * Uf - Kf^2, Uc' * Sc^2 * Uc - Kc^2
        W = Uf * V * Uc'
        S = [Sc zeros(nc, nf); zeros(nf, nc) Sf]
        R = [zeros(nc, nc) W'; W zeros(nf, nf)]
        Gp, Gm = [Uc; Uf * V] / sqrt(2), [Uc; -Uf * V] / sqrt(2)
        Pw = [Uc * Uc' zeros(nc, nf); zeros(nf, nc) Uf * Uf']
        anti, out, tot = graph_departure(Kf, Lf, Kc, Lc, V)
        @test anti ≈ opnorm(Gm' * S * Gp) rtol = 1e-8
        @test out ≈ opnorm((I - Pw) * S * Gp) rtol = 1e-7
        @test tot ≈ opnorm((I - Gp * Gp') * S * Gp) rtol = 1e-8
        @test window_defect(Kf, Lf, Kc, Lc, V) ≈ opnorm(Sf * W - W * Sc) rtol = 1e-8
        @test opnorm(S * R - R * S) ≈ opnorm(Sf * W - W * Sc) rtol = 1e-10
    end
end

end
