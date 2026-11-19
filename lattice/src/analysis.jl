"""
    mean_err(x; c = 6) -> (mean, error, τ_int)

Mean of a Monte Carlo time series with an error that accounts for
autocorrelations. τ_int is summed up to the first window W ≥ c τ_int(W)
(Madras–Sokal automatic windowing), and error = √(2 τ_int Γ(0) / N).
"""
function mean_err(x::AbstractVector{<:Real}; c::Real = 6)
    N = length(x)
    N >= 2 || throw(ArgumentError("need at least two measurements"))
    m = sum(x) / N
    d = x .- m
    Γ0 = sum(abs2, d) / N
    Γ0 == 0 && return (m, 0.0, 0.5)
    τ = 0.5
    for t in 1:(N ÷ 2)
        Γt = sum(d[i] * d[i+t] for i in 1:N-t) / (N - t)
        τ += Γt / Γ0
        t >= c * τ && break
    end
    τ = max(τ, 0.5)
    return (m, sqrt(2τ * Γ0 / N), τ)
end
