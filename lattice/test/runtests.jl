using Test, Random, LinearAlgebra
using LatticeYM
using LatticeYM: SUBGROUPS, lmul_su2, su2_near_identity, sample_a0

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

end
