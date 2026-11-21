# Glueball and static-potential measurements, and the correlator analysis used
# with them. Direction D is Euclidean time throughout.

"""
    slice_plaquettes(g) -> Vector

For each time slice t, the average of Re tr P_ij / 3 over the spatial
plaquettes (i < j < D) of that slice. Summed over a slice, the three
orientations form the A₁⁺⁺ representation of the cubic group, so this is a
zero-momentum 0⁺⁺ glueball operator. Equal to 1 on the unit configuration.
"""
function slice_plaquettes(g::GaugeField{D}) where {D}
    lat = g.lat
    T = lat.dims[D]
    O = zeros(T)
    for x in 1:nsites(lat)
        s = 0.0
        for i in 1:D-1, j in i+1:D-1
            s += retr(plaquette(g, i, j, x))
        end
        O[lat.coords[D, x]] += s
    end
    nplanes = (D - 1) * (D - 2) ÷ 2
    return O ./ (3 * nplanes * (nsites(lat) ÷ T))
end

"""
    glueball_operators(g, ρ, levels) -> Matrix

Rows: the 0⁺⁺ slice operator after `levels[k]` steps of spatial stout smearing
with parameter ρ (levels ascending). Columns: time slices.
"""
function glueball_operators(g::GaugeField{D}, ρ::Real, levels::AbstractVector{<:Integer}) where {D}
    issorted(levels) || throw(ArgumentError("smearing levels must be ascending"))
    out = zeros(length(levels), g.lat.dims[D])
    h, done = g, 0
    for (k, n) in enumerate(levels)
        h = stout_spatial(h, ρ, n - done)
        done = n
        out[k, :] = slice_plaquettes(h)
    end
    return out
end

"""
    wilson_loops_rt(g, Rmax, Tmax) -> Matrix

W[R, T] = average of Re tr / 3 of R×T loops with the R side along a spatial
direction (averaged over the D−1 of them) and the T side along time. Straight
lines are built once per site, so all loops cost little more than one. Agrees
with `wilson_loop(g, μ, D, R, T)` averaged over μ. Untwisted lattices only.
"""
function wilson_loops_rt(g::GaugeField{D}, Rmax::Int, Tmax::Int) where {D}
    lat, U = g.lat, g.U
    all(==(1), lat.z) || throw(ArgumentError("wilson_loops_rt supports untwisted lattices only"))
    (all(Rmax < lat.dims[μ] for μ in 1:D-1) && Tmax < lat.dims[D]) ||
        throw(ArgumentError("loops must fit inside the lattice"))
    V = nsites(lat)
    # line[n][x] = U(x) U(x+ê) ... (n links) along direction e
    function lines(e, nmax)
        L = [Vector{Mat3}(undef, V) for _ in 1:nmax]
        Threads.@threads for x in 1:V
            P = ONE3
            y = x
            for n in 1:nmax
                P = P * U[e, y]
                y = lat.fwd[e, y]
                L[n][x] = P
            end
        end
        return L
    end
    shift(x, e, n) = (y = x; for _ in 1:n; y = lat.fwd[e, y]; end; y)
    tl = lines(D, Tmax)
    W = zeros(Rmax, Tmax)
    for μ in 1:D-1
        sl = lines(μ, Rmax)
        chunks = collect(Iterators.partition(1:V, cld(V, 16)))   # fixed chunks: thread-count independent
        acc = zeros(Rmax, Tmax, length(chunks))
        Threads.@threads for c in eachindex(chunks)
            for x in chunks[c], R in 1:Rmax
                xR = shift(x, μ, R)
                for T in 1:Tmax
                    xT = shift(x, D, T)
                    # x → x+Rμ̂ → x+Rμ̂+Tt̂ → x+Tt̂ → x
                    acc[R, T, c] += retr(sl[R][x] * tl[T][xR], sl[R][xT]' * tl[T][x]')
                end
            end
        end
        W .+= dropdims(sum(acc; dims = 3); dims = 3)
    end
    return W ./ (3 * (D - 1) * V)
end

# ---------------------------------------------------------------------------
# Analysis

"""
    bin_means(x, b)

Means of consecutive blocks of `b` rows of `x` (along the first dimension),
dropping an incomplete final block. Binning reduces autocorrelation before a
jackknife.
"""
function bin_means(x::AbstractArray, b::Int)
    n = size(x, 1) ÷ b
    n >= 2 || throw(ArgumentError("need at least two bins"))
    return cat((sum(selectdim(x, 1, (i-1)*b+1:i*b); dims = 1) ./ b for i in 1:n)...; dims = 1)
end

"""
    jackknife(f, x) -> (estimate, error)

Delete-one jackknife over the rows (first dimension) of `x`. `f` maps an array
of rows to a number or an array; the estimate is `f(x)` and the error has the
usual (n−1)/n factor.
"""
function jackknife(f, x::AbstractArray)
    n = size(x, 1)
    full = f(x)
    samples = [f(selectdim(x, 1, [j for j in 1:n if j != i])) for i in 1:n]
    mean_s = sum(samples) / n
    var = sum(s -> (s .- mean_s) .^ 2, samples) .* ((n - 1) / n)
    return full, sqrt.(var)
end

"""
    correlator_matrix(O) -> C

O has size (nconf, nops, T). Returns C[k, l, t] = ⟨O_k(t₀+t) O_l(t₀)⟩ − ⟨O_k⟩⟨O_l⟩,
averaged over t₀ and over configurations, symmetrized in (k, l), for
t = 0, …, T−1 (index t+1).
"""
function correlator_matrix(O::AbstractArray{<:Real,3})
    nconf, nops, T = size(O)
    m = dropdims(sum(O; dims = (1, 3)); dims = (1, 3)) ./ (nconf * T)
    C = zeros(nops, nops, T)
    for c in 1:nconf, t in 0:T-1, t0 in 1:T
        t1 = mod1(t0 + t, T)
        for k in 1:nops, l in 1:nops
            C[k, l, t+1] += O[c, k, t1] * O[c, l, t0]
        end
    end
    C ./= nconf * T
    for t in 1:T
        C[:, :, t] .-= m * m'
        C[:, :, t] .= (C[:, :, t] + C[:, :, t]') ./ 2
    end
    return C
end

"""
    correlator_products(O) -> (P, M)

Per-configuration ingredients of `correlator_matrix`, so that they can be
binned before the vacuum subtraction: P[c, k, l, t+1] = (1/T) Σ_{t₀}
O[c, k, t₀+t] O[c, l, t₀] and M[c, k] = (1/T) Σ_t O[c, k, t]. Then
`correlator_from_means(mean P, mean M)` equals `correlator_matrix(O)`.
"""
function correlator_products(O::AbstractArray{<:Real,3})
    nconf, nops, T = size(O)
    P = zeros(nconf, nops, nops, T)
    M = dropdims(sum(O; dims = 3); dims = 3) ./ T
    for c in 1:nconf, t in 0:T-1, t0 in 1:T
        t1 = mod1(t0 + t, T)
        for k in 1:nops, l in 1:nops
            P[c, k, l, t+1] += O[c, k, t1] * O[c, l, t0]
        end
    end
    return P ./ T, M
end

"""
    correlator_from_means(P, M)

C[k, l, t] = P[k, l, t] − M[k] M[l], symmetrized in (k, l); `P` and `M` are
averages of the arrays from `correlator_products` over configurations.
"""
function correlator_from_means(P::AbstractArray{<:Real,3}, M::AbstractVector{<:Real})
    C = similar(P, Float64)
    for t in axes(P, 3)
        A = P[:, :, t] .- M * M'
        C[:, :, t] .= (A + A') ./ 2
    end
    return C
end

"""
    cosh_mass(r, t, T)

Solve r = cosh(m (t − T/2)) / cosh(m (t + 1 − T/2)) for m ≥ 0, the effective
mass of a periodic correlator from the ratio C(t)/C(t+1). Returns NaN when the
ratio admits no solution (r ≤ 1, noise).
"""
function cosh_mass(r::Real, t::Real, T::Real)
    (isfinite(r) && r > 1) || return NaN
    f(m) = cosh(m * (t - T / 2)) / cosh(m * (t + 1 - T / 2)) - r
    lo, hi = 0.0, 1.0
    while f(hi) < 0
        hi *= 2
        hi > 50 && return NaN
    end
    for _ in 1:200
        mid = (lo + hi) / 2
        f(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
    effective_masses(c) -> Vector

Effective masses m(t), t = 0, …, T/2 − 1, of a periodic correlator c (index t+1).
"""
effective_masses(c::AbstractVector) =
    [cosh_mass(c[t+1] / c[t+2], t, length(c)) for t in 0:length(c)÷2-1]

"""
    gevp(C, t0, t) -> (λ, V)

Generalized eigenproblem C(t) v = λ C(t₀) v for the symmetric correlator matrix
`C` from `correlator_matrix` (time index t+1). Eigenvalues are in descending
order; columns of V are normalized with vᵀ C(t₀) v = 1. Without periodic
images, λ_n = exp(−E_n (t − t₀)) when as many states as operators contribute.
"""
function gevp(C::AbstractArray{<:Real,3}, t0::Int, t::Int)
    A = Symmetric(C[:, :, t+1])
    B = Symmetric(C[:, :, t0+1])
    F = eigen(A, B)
    p = sortperm(F.values; rev = true)
    return F.values[p], F.vectors[:, p]
end
