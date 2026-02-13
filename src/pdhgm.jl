# Primal dual hybrid gradient method, based on:
# https://doi.org/10.1016/j.jmr.2017.05.010

using LinearAlgebra

function PDHGM(K::AbstractMatrix, s::AbstractVector, α::Real; tol=10, τ=10 , σ=0.1)

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

function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::pdhgm)

    f = PDHGM(K, g, α, τ=solver.τ, σ=solver.σ)

    r = K * f - g

    return f, r
end
