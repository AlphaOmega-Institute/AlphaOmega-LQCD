"""
    mean_plaquette(g)

Average of Re tr(z P) / 3 over all plaquettes. The action density is
⟨S⟩ = 1 − mean_plaquette(g).
"""
function mean_plaquette(g::GaugeField{D}) where {D}
    s = 0.0
    for x in 1:nsites(g.lat), μ in 1:D, ν in μ+1:D
        s += plaquette_retr(g, μ, ν, x)
    end
    return s / (3 * nsites(g.lat) * D * (D - 1) ÷ 2)
end

"""
    plane_plaquettes(g)

D×D symmetric matrix of plaquette averages per (μ, ν) plane; zero diagonal.
"""
function plane_plaquettes(g::GaugeField{D}) where {D}
    P = zeros(D, D)
    for x in 1:nsites(g.lat), μ in 1:D, ν in μ+1:D
        P[μ, ν] += plaquette_retr(g, μ, ν, x)
    end
    P ./= 3 * nsites(g.lat)
    return P + P'
end

"""
    corner_and_bulk(g, μ, ν) -> (corner, rest)

Average Re tr(z P) / 3 over the twisted corner plaquettes of the (μ, ν) plane
and over all other plaquettes of that plane. The twist can be moved around by a
change of variables, so the two agree in expectation for a correct twist.
"""
function corner_and_bulk(g::GaugeField, μ::Int, ν::Int)
    lat = g.lat
    sc, nc, sb, nb = 0.0, 0, 0.0, 0
    for x in 1:nsites(lat)
        v = plaquette_retr(g, μ, ν, x) / 3
        if is_corner(lat, μ, ν, x)
            sc += v; nc += 1
        else
            sb += v; nb += 1
        end
    end
    return sc / nc, sb / nb
end

"""
    polyakov_loop(g; dir = D)

Spatial average of tr(∏ U_dir) / 3 along closed lines in direction `dir`.
"""
function polyakov_loop(g::GaugeField{D}; dir::Int = D) where {D}
    lat = g.lat
    total = 0.0im
    n = 0
    for x in 1:nsites(lat)
        lat.coords[dir, x] == 1 || continue
        P = ONE3
        y = x
        for _ in 1:lat.dims[dir]
            P = P * g.U[dir, y]
            y = lat.fwd[dir, y]
        end
        total += tr(P)
        n += 1
    end
    return total / (3n)
end

"""
    wilson_loop(g, μ, ν, R, T)

Average of Re tr(W) / 3 over all R×T rectangular loops in the (μ, ν) plane,
with R < L_μ and T < L_ν. A loop that encloses the twisted corner plaquette
is multiplied by that plaquette's phase, so the 1×1 loop equals the plaquette.
"""
function wilson_loop(g::GaugeField{D}, μ::Int, ν::Int, R::Int, T::Int) where {D}
    lat = g.lat
    (R < lat.dims[μ] && T < lat.dims[ν]) || throw(ArgumentError("loop must fit inside the lattice"))
    total = 0.0
    for x in 1:nsites(lat)
        W = ONE3
        y = x
        for _ in 1:R
            W = W * g.U[μ, y]; y = lat.fwd[μ, y]
        end
        for _ in 1:T
            W = W * g.U[ν, y]; y = lat.fwd[ν, y]
        end
        for _ in 1:R
            y = lat.bwd[μ, y]; W = W * g.U[μ, y]'
        end
        for _ in 1:T
            y = lat.bwd[ν, y]; W = W * g.U[ν, y]'
        end
        encloses = mod(lat.dims[μ] - lat.coords[μ, x], lat.dims[μ]) < R &&
                   mod(lat.dims[ν] - lat.coords[ν, x], lat.dims[ν]) < T
        total += real((encloses ? lat.z[μ, ν] : one(ComplexF64)) * tr(W))
    end
    return total / (3 * nsites(lat))
end
