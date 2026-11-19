# 3×3 complex matrices stored column-major in an NTuple, so they live on the
# stack and the small products below compile to straight-line code.

struct Mat3
    e::NTuple{9,ComplexF64}
end

@inline Base.getindex(A::Mat3, i::Integer, j::Integer) = @inbounds A.e[i+3(j-1)]
@inline mat3(f::F) where {F} = Mat3(ntuple(k -> f((k - 1) % 3 + 1, (k - 1) ÷ 3 + 1), Val(9)))

Mat3(M::AbstractMatrix) = mat3((i, j) -> ComplexF64(M[i, j]))
Base.Matrix(A::Mat3) = [A[i, j] for i in 1:3, j in 1:3]

const ONE3 = mat3((i, j) -> i == j ? 1.0 + 0.0im : 0.0 + 0.0im)
const ZERO3 = Mat3(ntuple(_ -> 0.0 + 0.0im, Val(9)))
Base.one(::Type{Mat3}) = ONE3
Base.zero(::Type{Mat3}) = ZERO3

@inline Base.:+(A::Mat3, B::Mat3) = Mat3(map(+, A.e, B.e))
@inline Base.:-(A::Mat3, B::Mat3) = Mat3(map(-, A.e, B.e))
@inline Base.:*(s::Number, A::Mat3) = Mat3(map(x -> s * x, A.e))
@inline Base.:*(A::Mat3, B::Mat3) =
    mat3((i, j) -> A[i, 1] * B[1, j] + A[i, 2] * B[2, j] + A[i, 3] * B[3, j])
@inline Base.adjoint(A::Mat3) = mat3((i, j) -> conj(A[j, i]))
@inline LinearAlgebra.tr(A::Mat3) = A[1, 1] + A[2, 2] + A[3, 3]
Base.isapprox(A::Mat3, B::Mat3; kwargs...) = isapprox(collect(A.e), collect(B.e); kwargs...)

@inline retr(A::Mat3) = real(tr(A))

"""
    retr(A, B)

Re tr(A B), without forming the product.
"""
@inline function retr(A::Mat3, B::Mat3)
    s = 0.0
    @inbounds for i in 1:3, k in 1:3
        s += real(A[i, k] * B[k, i])
    end
    return s
end

function LinearAlgebra.det(A::Mat3)
    return A[1, 1] * (A[2, 2] * A[3, 3] - A[2, 3] * A[3, 2]) -
           A[1, 2] * (A[2, 1] * A[3, 3] - A[2, 3] * A[3, 1]) +
           A[1, 3] * (A[2, 1] * A[3, 2] - A[2, 2] * A[3, 1])
end

"""
    su3_deviation(A)

Distance of `A` from SU(3): the larger of max|A A† − 1| and |det A − 1|.
"""
su3_deviation(A::Mat3) = max(maximum(abs, (A * A' - ONE3).e), abs(det(A) - 1))

"""
    reunitarize(A)

Project `A` onto SU(3): Gram–Schmidt on the first two rows, then the third row
is fixed as conj(r₁ × r₂), which makes det = 1 exactly.
"""
@inline function reunitarize(A::Mat3)
    # Each row gets a fresh name: reassigning a variable that the closure below
    # captures would box it and allocate.
    a = (A[1, 1], A[1, 2], A[1, 3])
    r1 = a ./ sqrt(sum(abs2, a))
    b = (A[2, 1], A[2, 2], A[2, 3])
    c = conj(r1[1]) * b[1] + conj(r1[2]) * b[2] + conj(r1[3]) * b[3]
    b⊥ = b .- c .* r1
    r2 = b⊥ ./ sqrt(sum(abs2, b⊥))
    r3 = conj.((r1[2] * r2[3] - r1[3] * r2[2],
                r1[3] * r2[1] - r1[1] * r2[3],
                r1[1] * r2[2] - r1[2] * r2[1]))
    return mat3((i, j) -> i == 1 ? r1[j] : i == 2 ? r2[j] : r3[j])
end

"""
    haar_su3(rng)

Haar-random SU(3) matrix. Gram–Schmidt on Gaussian rows is right-invariant,
so the result is Haar distributed.
"""
haar_su3(rng::AbstractRNG) = reunitarize(Mat3(ntuple(_ -> randn(rng, ComplexF64), Val(9))))

# ---------------------------------------------------------------------------
# SU(2) subgroups. An SU(2) element [a b; -conj(b) conj(a)] with |a|²+|b|²=1
# is stored as the pair (a, b); in terms of a₀ + i a⃗·σ⃗ this is
# a = a₀ + i a₃, b = a₂ + i a₁.

const SUBGROUPS = ((1, 2), (1, 3), (2, 3))

@inline su2mul(x::NTuple{2,ComplexF64}, y::NTuple{2,ComplexF64}) =
    (x[1] * y[1] - x[2] * conj(y[2]), x[1] * y[2] + x[2] * conj(y[1]))

"""
    lmul_su2(a, b, i, j, U)

R U, where R is the SU(2) element (a, b) embedded in rows/columns (i, j) of
the identity. Only rows i and j of U change.
"""
@inline function lmul_su2(a::ComplexF64, b::ComplexF64, i::Int, j::Int, U::Mat3)
    return mat3((r, c) -> r == i ? a * U[i, c] + b * U[j, c] :
                          r == j ? -conj(b) * U[i, c] + conj(a) * U[j, c] : U[r, c])
end

"""
    subgroup_projection(W, i, j) -> (p, q, k)

For R in the (i, j) subgroup, Re tr(R W) = Re(a p) + Re(b q) + const, with
p = W_ii + conj(W_jj) and q = W_ji − conj(W_ij). Writing this as
k Re tr(r v) with v ∈ SU(2) gives k = √(|p|² + |q|²)/2 and
v† = (conj(p), conj(q)) / 2k.
"""
@inline function subgroup_projection(W::Mat3, i::Int, j::Int)
    p = W[i, i] + conj(W[j, j])
    q = W[j, i] - conj(W[i, j])
    return p, q, sqrt(abs2(p) + abs2(q)) / 2
end

"""Random SU(2) element with a₀ fixed and a⃗ uniform on the sphere of radius √(1−a₀²)."""
@inline function su2_with_a0(a0::Float64, rng::AbstractRNG)
    r = sqrt(max(0.0, 1 - a0^2))
    cosθ = 2rand(rng) - 1
    sinθ = sqrt(max(0.0, 1 - cosθ^2))
    φ = 2π * rand(rng)
    a1, a2, a3 = r * sinθ * cos(φ), r * sinθ * sin(φ), r * cosθ
    return (complex(a0, a3), complex(a2, a1))
end

"""Random SU(2) element at distance ε from the identity; symmetric under inversion."""
@inline function su2_near_identity(ε::Float64, rng::AbstractRNG)
    n1, n2, n3 = randn(rng), randn(rng), randn(rng)
    s = ε / sqrt(n1^2 + n2^2 + n3^2)
    return (complex(sqrt(1 - ε^2), s * n3), complex(s * n2, s * n1))
end

"""
    sample_a0(α, rng)

Sample a₀ ∈ [−1, 1] with density ∝ √(1−a₀²) exp(α a₀). Uses Kennedy–Pendleton
for α > 2 and rejection from the semicircle otherwise.
"""
function sample_a0(α::Float64, rng::AbstractRNG)
    if α > 2
        while true
            r1, r2, r3, r4 = 1 - rand(rng), rand(rng), 1 - rand(rng), rand(rng)
            λ2 = -(log(r1) + cospi(2r2)^2 * log(r3)) / (2α)
            r4^2 <= 1 - λ2 && return 1 - 2λ2
        end
    else
        while true
            # The x coordinate of a uniform point in the unit disk has density ∝ √(1−x²).
            x, y = 2rand(rng) - 1, 2rand(rng) - 1
            x^2 + y^2 > 1 && continue
            rand(rng) <= exp(α * (x - 1)) && return x
        end
    end
end
