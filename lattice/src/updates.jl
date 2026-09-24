abstract type LinkUpdate end

"""
    Metropolis(; ε = 0.3, nhit = 4)

Metropolis with `nhit` hits per link. Each proposal multiplies the link by a
random SU(2)-subgroup element at distance ε from the identity; the proposal
distribution is symmetric under inversion, so detailed balance holds.
"""
struct Metropolis <: LinkUpdate
    ε::Float64
    nhit::Int
end
Metropolis(; ε = 0.3, nhit = 4) = Metropolis(ε, nhit)

"""
    Heatbath()

Cabibbo–Marinari heatbath: one exact SU(2) heatbath step in each of the three
SU(2) subgroups, using Kennedy–Pendleton sampling.
"""
struct Heatbath <: LinkUpdate end

"""
    Overrelaxation()

Microcanonical overrelaxation in the three SU(2) subgroups. It leaves the
action unchanged, so it is not ergodic on its own; alternate it with `Heatbath`.
"""
struct Overrelaxation <: LinkUpdate end

hits(alg::Metropolis) = alg.nhit
hits(::LinkUpdate) = 1

# Each (μ, parity) batch is split into this many chunks, each with its own RNG
# seeded from the caller's RNG, so results do not depend on the thread count.
const NCHUNKS = 16

"""
    sweep!(g, β, alg, rng) -> acceptance rate

Update every link once with `alg`. Links are processed in batches of fixed
direction μ and site parity: the staple of U_μ(x) only contains μ-links at
x ± ν̂, which have the other parity, so a batch can be updated in parallel.
"""
function sweep!(g::GaugeField{D}, β::Real, alg::LinkUpdate, rng::AbstractRNG) where {D}
    lat = g.lat
    βf = Float64(β)
    accepted = 0
    tried = 0
    for μ in 1:D, p in 1:2
        sites = lat.parity_sites[p]
        n = length(sites)
        seeds = rand(rng, UInt64, NCHUNKS)
        counts = zeros(Int, NCHUNKS)
        run_chunk = c -> begin
            lrng = Xoshiro(seeds[c])
            acc = 0
            @inbounds for i in (div((c - 1) * n, NCHUNKS) + 1):div(c * n, NCHUNKS)
                x = sites[i]
                Unew, a = update_link(alg, g.U[μ, x], staple(g, μ, x), βf, lrng)
                g.U[μ, x] = Unew
                acc += a
            end
            counts[c] = acc
        end
        if Threads.nthreads() == 1 || n < 1024
            foreach(run_chunk, 1:NCHUNKS)
        else
            @sync for c in 1:NCHUNKS
                Threads.@spawn run_chunk(c)
            end
        end
        accepted += sum(counts)
        tried += n
    end
    return accepted / (tried * hits(alg))
end

# The local weight of U is exp((β/3) Re tr(U A)).

function update_link(alg::Metropolis, U::Mat3, A::Mat3, β::Float64, rng::AbstractRNG)
    s_old = retr(U, A)
    acc = 0
    for _ in 1:alg.nhit
        i, j = SUBGROUPS[rand(rng, 1:3)]
        a, b = su2_near_identity(alg.ε, rng)
        Unew = lmul_su2(a, b, i, j, U)
        s_new = retr(Unew, A)
        if rand(rng) < exp(β / 3 * (s_new - s_old))
            U, s_old = Unew, s_new
            acc += 1
        end
    end
    return reunitarize(U), acc
end

function update_link(::Heatbath, U::Mat3, A::Mat3, β::Float64, rng::AbstractRNG)
    for (i, j) in SUBGROUPS
        # Re tr(R U A) = k Re tr(r v) + const; draw s = r v from exp((2βk/3) s₀).
        p, q, k = subgroup_projection(U * A, i, j)
        s = su2_with_a0(sample_a0(2β * k / 3, rng), rng)
        a, b = k > 0 ? su2mul(s, (conj(p) / 2k, conj(q) / 2k)) : s
        U = lmul_su2(a, b, i, j, U)
    end
    return reunitarize(U), 1
end

function update_link(::Overrelaxation, U::Mat3, A::Mat3, β::Float64, rng::AbstractRNG)
    for (i, j) in SUBGROUPS
        # Map r v = v to r v = v†: Re tr is unchanged and the map is an involution.
        p, q, k = subgroup_projection(U * A, i, j)
        k > 0 || continue
        vdag = (conj(p) / 2k, conj(q) / 2k)
        a, b = su2mul(vdag, vdag)
        U = lmul_su2(a, b, i, j, U)
    end
    return reunitarize(U), 1
end
