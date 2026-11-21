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

"""
    su3_project(A)

Project a 3×3 matrix onto SU(3) through its polar decomposition:
W = A (A†A)^{-1/2}, then W times a cube root of conj(det W). For G, H in SU(3),
`su3_project(G A H') = G su3_project(A) H'`, so the projection is covariant at
both ends of a path. (`reunitarize`, row Gram–Schmidt, is covariant only under
right multiplication.) `A` must be invertible.
"""
function su3_project(A::Mat3)
    M = Matrix(A)
    F = eigen(Hermitian(M' * M))
    W = M * F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors'
    return Mat3(cis(-angle(det(W)) / 3) * W)
end

"""
    block_spatial(g, α) -> GaugeField

One blocking step in the spatial directions. The coarse lattice has spatial
extents L_i/2 (which must be even) and the same time extent; coarse site X sits
on fine site x with x_i = 2X_i − 1. A coarse spatial link is the SU(3)
projection of the straight path U_μ(x) U_μ(x+μ̂) plus α times the 2(D−2) bent
paths of length four from x to x+2μ̂ through x ± ν̂ (ν spatial, ν ≠ μ). Every path
runs from x to x+2μ̂, so with `su3_project` the blocked links transform under
the gauge transformation restricted to coarse sites. Time links are copied
from the fine site and are not blocked. Untwisted lattices only.
"""
function block_spatial(g::GaugeField{D}, α::Real) where {D}
    lat, U = g.lat, g.U
    all(==(1), lat.z) || throw(ArgumentError("block_spatial supports untwisted lattices only"))
    cdims = ntuple(d -> d == D ? lat.dims[d] : lat.dims[d] ÷ 2, D)
    all(d -> iseven(lat.dims[d]) && iseven(cdims[d]), 1:D-1) ||
        throw(ArgumentError("spatial extents must be multiples of 4"))
    clat = Lattice(cdims)
    li = LinearIndices(lat.dims)
    Uc = Matrix{Mat3}(undef, D, nsites(clat))
    αf = Float64(α)
    Threads.@threads for X in 1:nsites(clat)
        x = li[CartesianIndex(ntuple(d -> d == D ? clat.coords[d, X] : 2clat.coords[d, X] - 1, D))]
        Uc[D, X] = U[D, x]
        for μ in 1:D-1
            xμ = lat.fwd[μ, x]
            x2μ = lat.fwd[μ, xμ]
            A = U[μ, x] * U[μ, xμ]
            for ν in 1:D-1
                ν == μ && continue
                xp = lat.fwd[ν, x]                     # through x + ν̂
                A = A + αf * (U[ν, x] * U[μ, xp] * U[μ, lat.fwd[μ, xp]] * U[ν, x2μ]')
                xm = lat.bwd[ν, x]                     # through x − ν̂
                A = A + αf * (U[ν, xm]' * U[μ, xm] * U[μ, lat.fwd[μ, xm]] * U[ν, lat.bwd[ν, x2μ]])
            end
            Uc[μ, X] = su3_project(A)
        end
    end
    return GaugeField(clat, Uc)
end
