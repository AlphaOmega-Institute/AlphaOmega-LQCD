"""
    expm3(X)

Exponential of a 3×3 matrix by scaling and squaring with a Taylor polynomial of
order 12. For anti-Hermitian traceless `X` the result is in SU(3) to rounding.
The scaling depends only on the Frobenius norm of `X`, so `expm3(G X G†) =
G expm3(X) G†` for unitary `G`, up to rounding.
"""
@inline function expm3(X::Mat3)
    nrm = sqrt(sum(abs2, X.e))
    s = 0
    while nrm > 0.25
        nrm /= 2
        s += 1
    end
    Y = (1 / 2^s) * X
    E = ONE3
    term = ONE3
    for k in 1:12
        term = (1 / k) * (term * Y)
        E = E + term
    end
    for _ in 1:s
        E = E * E
    end
    return E
end

"""
    stout_spatial(g, ρ) -> GaugeField

One step of stout smearing (Morningstar–Peardon) applied to the spatial links,
directions 1 to D−1; the links in direction D (time) are copied unchanged. With
C_μ(x) = ρ Σ_ν [U_ν(x) U_μ(x+ν̂) U_ν(x+μ̂)† + U_ν(x−ν̂)† U_μ(x−ν̂) U_ν(x−ν̂+μ̂)],
summed over spatial ν ≠ μ, and Ω = C_μ(x) U_μ(x)†, the new link is
exp(X) U_μ(x), where X = (Ω − Ω†)/2 − tr(Ω − Ω†)/6 is the traceless
anti-Hermitian part. Ω transforms as G(x) Ω G(x)†, so smearing commutes with
gauge transformations exactly. Only untwisted lattices are supported.
"""
function stout_spatial(g::GaugeField{D}, ρ::Real) where {D}
    lat, U = g.lat, g.U
    all(==(1), lat.z) || throw(ArgumentError("stout_spatial supports untwisted lattices only"))
    D >= 3 || throw(ArgumentError("need at least two spatial directions"))
    Unew = copy(U)
    ρf = Float64(ρ)
    Threads.@threads for x in 1:nsites(lat)
        @inbounds for μ in 1:D-1
            C = ZERO3
            xpμ = lat.fwd[μ, x]
            for ν in 1:D-1
                ν == μ && continue
                xpν = lat.fwd[ν, x]
                xmν = lat.bwd[ν, x]
                xpμmν = lat.fwd[μ, xmν]
                C = C + U[ν, x] * U[μ, xpν] * U[ν, xpμ]' + U[ν, xmν]' * U[μ, xmν] * U[ν, xpμmν]
            end
            Ω = ρf * (C * U[μ, x]')
            A = Ω - Ω'
            X = 0.5 * A - (tr(A) / 6) * ONE3
            Unew[μ, x] = expm3(X) * U[μ, x]
        end
    end
    return GaugeField(lat, Unew)
end

"""
    stout_spatial(g, ρ, n)

`n` steps of spatial stout smearing.
"""
function stout_spatial(g::GaugeField, ρ::Real, n::Integer)
    h = g
    for _ in 1:n
        h = stout_spatial(h, ρ)
    end
    return h
end
