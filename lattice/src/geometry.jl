"""
    Lattice(dims; twist = ())

Periodic hypercubic lattice with extents `dims`, which must all be even so the
lattice splits into two checkerboard parities.

`twist` lists 't Hooft twists as `(μ, ν) => n` pairs, e.g.
`((1, 2) => 1, (3, 4) => 1)`. In each (μ, ν) plane-slice exactly one
plaquette, the one at x_μ = L_μ, x_ν = L_ν, is multiplied by the center
element z = exp(2πi n / 3); the reversed orientation carries conj(z).
"""
struct Lattice{D}
    dims::NTuple{D,Int}
    coords::Matrix{Int}               # coords[μ, x], 1-based
    fwd::Matrix{Int}                  # fwd[μ, x] = x + μ̂
    bwd::Matrix{Int}                  # bwd[μ, x] = x − μ̂
    parity_sites::NTuple{2,Vector{Int}}
    z::Matrix{ComplexF64}             # z[μ, ν] = conj(z[ν, μ])
end

function Lattice(dims::NTuple{D,Int}; twist = ()) where {D}
    all(L -> L >= 2 && iseven(L), dims) ||
        throw(ArgumentError("lattice extents must be even and ≥ 2, got $dims"))
    V = prod(dims)
    li = LinearIndices(dims)
    coords = Matrix{Int}(undef, D, V)
    fwd = similar(coords)
    bwd = similar(coords)
    for (x, I) in enumerate(CartesianIndices(dims))
        for μ in 1:D
            coords[μ, x] = I[μ]
            fwd[μ, x] = li[CartesianIndex(ntuple(d -> d == μ ? mod1(I[d] + 1, dims[d]) : I[d], D))]
            bwd[μ, x] = li[CartesianIndex(ntuple(d -> d == μ ? mod1(I[d] - 1, dims[d]) : I[d], D))]
        end
    end
    even = [x for x in 1:V if iseven(sum(@view coords[:, x]))]
    odd = [x for x in 1:V if isodd(sum(@view coords[:, x]))]

    z = ones(ComplexF64, D, D)
    for ((μ, ν), n) in twist
        (1 <= μ <= D && 1 <= ν <= D && μ != ν) ||
            throw(ArgumentError("invalid twist plane ($μ, $ν) for a $D-dimensional lattice"))
        z[μ, ν] = cispi(2n / 3)
        z[ν, μ] = conj(z[μ, ν])
    end
    return Lattice{D}(dims, coords, fwd, bwd, (even, odd), z)
end

nsites(lat::Lattice) = size(lat.coords, 2)
Base.ndims(::Lattice{D}) where {D} = D

"""Whether x is the twisted corner of its (μ, ν) plane-slice."""
@inline is_corner(lat::Lattice, μ::Int, ν::Int, x::Int) =
    @inbounds lat.coords[μ, x] == lat.dims[μ] && lat.coords[ν, x] == lat.dims[ν]

"""Twist phase carried by the plaquette P_μν(x)."""
@inline zfac(lat::Lattice, μ::Int, ν::Int, x::Int) =
    is_corner(lat, μ, ν, x) ? @inbounds(lat.z[μ, ν]) : one(ComplexF64)
