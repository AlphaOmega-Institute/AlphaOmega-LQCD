"""
    single_plaquette_average(β; n = 256)

Exact ⟨Re tr U / 3⟩ for a single SU(3) matrix with weight exp((β/3) Re tr U).

This is the exact infinite-volume plaquette of two-dimensional SU(3) lattice
gauge theory, and the leading term of the four-dimensional strong-coupling
expansion (corrections start at order u⁵). Computed with the Weyl integration
formula over the eigenphases θ₁, θ₂, θ₃ = −θ₁ − θ₂; the integrand is smooth
and periodic, so the trapezoid rule converges exponentially in `n`.
"""
function single_plaquette_average(β::Real; n::Int = 256)
    num = 0.0
    den = 0.0
    for i in 0:n-1, j in 0:n-1
        θ1, θ2 = 2π * i / n, 2π * j / n
        θ3 = -θ1 - θ2
        vandermonde = (2sin((θ1 - θ2) / 2))^2 * (2sin((θ1 - θ3) / 2))^2 * (2sin((θ2 - θ3) / 2))^2
        t = cos(θ1) + cos(θ2) + cos(θ3)
        w = vandermonde * exp(β / 3 * (t - 3))
        num += w * t / 3
        den += w
    end
    return num / den
end
