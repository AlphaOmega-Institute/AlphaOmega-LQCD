"""
    GaugeField(lat, U)

SU(3) link variables on `lat`; `U[μ, x]` is the link from x to x + μ̂.
"""
struct GaugeField{D}
    lat::Lattice{D}
    U::Matrix{Mat3}
end

cold_start(lat::Lattice) = GaugeField(lat, fill(ONE3, ndims(lat), nsites(lat)))
hot_start(lat::Lattice, rng::AbstractRNG) =
    GaugeField(lat, [haar_su3(rng) for _ in 1:ndims(lat), _ in 1:nsites(lat)])
Base.copy(g::GaugeField) = GaugeField(g.lat, copy(g.U))

"""
    plaquette(g, μ, ν, x)

P_μν(x) = U_μ(x) U_ν(x+μ̂) U_μ(x+ν̂)† U_ν(x)†, without the twist phase.
"""
@inline function plaquette(g::GaugeField, μ::Int, ν::Int, x::Int)
    lat = g.lat
    @inbounds return g.U[μ, x] * g.U[ν, lat.fwd[μ, x]] * g.U[μ, lat.fwd[ν, x]]' * g.U[ν, x]'
end

"""
    plaquette_retr(g, μ, ν, x)

Re tr(z_μν(x) P_μν(x)), the quantity that enters the action.
"""
@inline plaquette_retr(g::GaugeField, μ::Int, ν::Int, x::Int) =
    real(zfac(g.lat, μ, ν, x) * tr(plaquette(g, μ, ν, x)))

"""
    staple(g, μ, x)

Sum of staples A around U_μ(x), defined so that the plaquettes containing that
link contribute Σ Re tr(z P) = Re tr(U_μ(x) A). Twist phases are included;
`test/runtests.jl` checks this against a full recomputation of the action.
"""
@inline function staple(g::GaugeField{D}, μ::Int, x::Int) where {D}
    lat, U = g.lat, g.U
    A = ZERO3
    @inbounds begin
        xpμ = lat.fwd[μ, x]
        for ν in 1:D
            ν == μ && continue
            xpν = lat.fwd[ν, x]
            xmν = lat.bwd[ν, x]
            xpμmν = lat.fwd[μ, xmν]
            up = U[ν, xpμ] * U[μ, xpν]' * U[ν, x]'            # closes P_μν(x)
            dn = U[ν, xpμmν]' * U[μ, xmν]' * U[ν, xmν]        # closes P_μν(x − ν̂)†
            A = A + zfac(lat, μ, ν, x) * up + conj(zfac(lat, μ, ν, xmν)) * dn
        end
    end
    return A
end

"""
    action(g, β)

Wilson action S = β Σ_{x, μ<ν} (1 − Re tr(z P_μν(x)) / 3).
"""
function action(g::GaugeField{D}, β::Real) where {D}
    s = 0.0
    for x in 1:nsites(g.lat), μ in 1:D, ν in μ+1:D
        s += 1 - plaquette_retr(g, μ, ν, x) / 3
    end
    return β * s
end

"""
    gauge_transform!(g, G)

Apply U_μ(x) → G(x) U_μ(x) G(x+μ̂)†. Twist phases are central, so every
observable in this package is invariant.
"""
function gauge_transform!(g::GaugeField{D}, G::AbstractVector{Mat3}) where {D}
    lat = g.lat
    for x in 1:nsites(lat), μ in 1:D
        g.U[μ, x] = G[x] * g.U[μ, x] * G[lat.fwd[μ, x]]'
    end
    return g
end
