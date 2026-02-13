# Primal dual hybrid gradient method, based on:
# https://doi.org/10.1016/j.jmr.2017.05.010

using LinearAlgebra

function PDHGM(K::AbstractMatrix, s::AbstractVector, α::Real; tol=1e-5, τ=0.1 , σ=10)

    B = inv(LinearAlgebra.I + τ * α * K' * K)
    Y = zeros(size(K, 2))
    Ỹ = copy(Y)
    f = ones(size(K, 2))
    f̃ = copy(f)
    f_prev = deepcopy(f)
    ε = tol + 1

    while ε > tol
        Ỹ .= Y + σ * f̃
        Y .= Ỹ ./ max.(1, abs.(Ỹ))
        f .= B * (f̃ - τ .* Y + τ * α * K' * s)
        f .= max.(0, f)
        f̃ .= 2 .* f .- f_prev
        ε = norm(f - f_prev) / norm(f_prev)
        f_prev .= f
    end
    return f
end

function PDHGM2(K::AbstractMatrix, s::AbstractVector, α::Real; tol=1e-5, τ=0.1, σ=10)

    n = size(K, 2)

    Y, Ỹ, temp_vec = zeros.((n, n, n))
    f, f̃, f_prev  = ones.((n, n, n))
    
    Ks = (τ * α) .* (K' * s)   

    B = inv(I + τ * α * K' * K)

    ε = tol + 1.0

    while ε > tol

        @. Ỹ = Y + σ * f̃
        @. Y = Ỹ / max(1, abs(Ỹ))
        
        @. temp_vec = f̃ - (τ * Y) + Ks
        mul!(f, B, temp_vec) 
        @. f = max(0, f)
        
        num = 0.0
        den = 0.0
        @inbounds for i in eachindex(f)
            diff = f[i] - f_prev[i]
            num += diff * diff
            den += f_prev[i] * f_prev[i]
        end
        ε = sqrt(num) / sqrt(den)
        
        @. f̃ = 2 * f - f_prev
        copyto!(f_prev, f)
    end
    
    return f
end

function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::pdhgm)

    f = PDHGM2(K, g, α, tol=solver.tol, τ=solver.τ, σ=solver.σ)

    r = K * f - g

    return f, r
end
