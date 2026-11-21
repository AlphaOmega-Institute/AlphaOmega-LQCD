# Checks of the inequalities in spectral_transfer.tex on random finite-dimensional
# instances. Run: julia docs/notes/check_spectral_transfer.jl   (standard library only)

using LinearAlgebra, Random, Test

herm(rng, n) = (A = randn(rng, ComplexF64, n, n); (A + A') / 2)
isometry(rng, n, m) = Matrix(qr(randn(rng, ComplexF64, n, m)).Q)[:, 1:m]
distance(x, S) = minimum(abs.(x .- S))

"""Fine Hamiltonian built around a coarse one: H_f = W H_c W' + (1-P) B (1-P) + coupling."""
function instance(rng; n = 40, m = 8, Λ0 = 6.0, κ = 0.3)
    W = isometry(rng, n, m)
    Hc = Diagonal(sort([0.0; 1.0 .+ 3rand(rng, m - 1)]))           # simple ground at 0, gap ≥ 1
    Hc = Matrix(Hc)
    P = W * W'
    B = Λ0 * I + 2 * (x -> x' * x)(randn(rng, ComplexF64, n, n)) / n
    X = κ * herm(rng, n)
    Hf = W * Hc * W' + (I - P) * B * (I - P) + P * X * (I - P) + (I - P) * X * P
    return Hermitian(Hf), Hermitian(Hc), W
end

@testset "spectral transfer" begin
    rng = Xoshiro(2026)
    for trial in 1:200
        κ = rand(rng, (0.0, 0.05, 0.2, 0.5, 1.0))
        Hf, Hc, W = instance(rng; κ)
        P = W * W'
        δ = opnorm(Hf * W - W * Hc)
        ef, Vf = eigen(Hf)
        ec = eigvals(Hc)
        # Proposition 1(i): coarse levels shadow fine levels
        @test all(μ -> distance(μ, ef) <= δ + 1e-10, ec)
        # Proposition 1(ii): below the UV bound, fine eigenvectors lie near Ran W and shadow coarse levels
        Λ = eigmin(Hermitian((I - P) * Hf * (I - P) + 1e6 * P))    # compression to Ran(1-P)
        for (λ, φ) in zip(ef, eachcol(Vf))
            λ < Λ - 1e-8 || continue
            η = δ / (Λ - λ)
            @test norm((I - P) * φ) <= η * norm(P * φ) + 1e-10
            @test distance(λ, ec) <= δ * sqrt(1 + η^2) + 1e-10
        end
        # Proposition 2: gap transfer
        Ec, Δc = ec[1], ec[2] - ec[1]
        λ0, λ1 = ef[1], ef[2]
        if λ1 < Λ - 1e-8 && Ec + Δc < Λ
            η1 = δ / (Λ - λ1)
            ε1 = δ * sqrt(1 + η1^2)
            κ2 = sqrt(2) * (δ + ε1) / sqrt(max(1 - 2η1^2, 1e-12))
            if 2η1^2 < 1 && κ2 < Δc
                @test λ0 <= Ec + δ + 1e-10
                @test λ1 - λ0 >= Δc - δ - ε1 - 1e-10
            end
        end
    end
    # the restricted version: a partial isometry W Q_E with the defect on the window only
    rng = Xoshiro(7)
    for trial in 1:100
        Hf, Hc, W = instance(rng; κ = 0.3)
        ec, Vc = eigen(Hc)
        QE = Vc[:, 1:3] * Vc[:, 1:3]'                     # lowest three coarse levels
        WE = W * QE
        δE = opnorm(Hf * WE - WE * Hc)
        @test δE <= opnorm(Hf * W - W * Hc) + 1e-10
        PE = WE * WE'
        ΛE = eigmin(Hermitian((I - PE) * Hf * (I - PE) + 1e6 * PE))
        ef, Vf = eigen(Hf)
        for (λ, φ) in zip(ef, eachcol(Vf))
            λ < ΛE - 1e-8 || continue
            η = δE / (ΛE - λ)
            @test norm((I - PE) * φ) <= η * norm(PE * φ) + 1e-10
            @test distance(λ, ec[1:3]) <= δE * sqrt(1 + η^2) + 1e-10
        end
    end
    # sharp residual bound of Remark 5: Var_ψ(H) ≥ g² f (1 − f)
    rng = Xoshiro(11)
    for _ in 1:200
        n = 8
        E = sort([0.0; 1.0 .+ 3rand(rng, n - 1)])
        g = E[2] - E[1]
        ψ = randn(rng, ComplexF64, n); ψ /= norm(ψ)
        f = abs2(ψ[1])
        H = Diagonal(E)
        δ = norm((H - real(ψ' * H * ψ) * I) * ψ)
        @test δ >= g * sqrt(f * (1 - f)) - 1e-12
    end
    for f in (0.1, 0.5, 0.93)                                   # equality for two levels
        ψ = [sqrt(f), sqrt(1 - f)]; H = Diagonal([0.0, 1.7])
        @test norm((H - (ψ' * H * ψ) * I) * ψ) ≈ 1.7 * sqrt(f * (1 - f))
    end
    # Proposition 4: two operator windows, with the alignment induced by an operator correspondence
    rng = Xoshiro(12)
    psd(rng, n) = (A = randn(rng, ComplexF64, n, n); M = A * A'; Hermitian(M / opnorm(M)))
    isqrt(G) = (F = eigen(Hermitian(G)); F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors')
    msqrt(G) = (F = eigen(Hermitian(G)); F.vectors * Diagonal(sqrt.(F.values)) * F.vectors')
    polar_u(A) = (F = svd(A); F.U * F.Vt)
    for trial in 1:300
        nf, nc, k = 30, 12, rand(rng, 1:4)
        Sf, Sc = psd(rng, nf), psd(rng, nc)
        F = randn(rng, ComplexF64, nf, k)                         # fine operator states f_j
        Cs = randn(rng, ComplexF64, nc, k)                        # coarse operator states c_j
        Gf, Gc = F' * F, Cs' * Cs
        Uf, Uc = F * isqrt(Gf), Cs * isqrt(Gc)
        V = polar_u(msqrt(Gf) * isqrt(Gc))
        W = Uf * V * Uc'
        # W is the polar part of the operator-state correspondence W0 c_j = f_j
        W0 = Uf * msqrt(Gf) * isqrt(Gc) * Uc'
        @test W0 * Cs ≈ F
        @test W ≈ polar_u(W0 * Uc) * Uc' atol = 1e-8
        # K and L from correlators C(t) = states' S^t states
        Cf(t) = F' * Sf^t * F
        Cc(t) = Cs' * Sc^t * Cs
        Kf = isqrt(Gf) * Cf(1) * isqrt(Gf); Lf = isqrt(Gf) * Cf(2) * isqrt(Gf) - Kf^2
        Kc = isqrt(Gc) * Cc(1) * isqrt(Gc); Lc = isqrt(Gc) * Cc(2) * isqrt(Gc) - Kc^2
        @test Kf ≈ Uf' * Sf * Uf atol = 1e-10
        @test Lf ≈ Uf' * Sf * (I - Uf * Uf') * Sf * Uf atol = 1e-10
        @test eigmin(Hermitian(Lc)) >= -1e-10
        lhs = opnorm(Sf * W - W * Sc)
        rhs = opnorm(Kf * V - V * Kc) + max(sqrt(opnorm(Lf)), sqrt(opnorm(Lc)))
        @test lhs <= rhs + 1e-10
    end
end
