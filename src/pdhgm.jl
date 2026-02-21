# Primal dual hybrid gradient method, based on:
# https://doi.org/10.1016/j.jmr.2017.05.010

using LinearAlgebra

export pdhgm
"""
    pdhgm(σ, τ, tol)
Primal dual hybrid gradient method for L1 regularization, 
following [this paper](https://doi.org/10.1016/j.jmr.2017.05.010)
[Reci2017](@cite).

It can be used as a "solver" for the invert function.

Positional (keyowrd) arguments:

- `sigma` (default value `10`)
- `tau` (default value `0.1`)
- `tol` (default value `1e-5`)

The particular choice of σ and τ is heuristic. 
The parameters σ and τ are step size parameters, which control
convergence and stability of the algorithm. Convergence is guaranteed 
when `τσ ≤ 1`.

The best values of σ and τ will depend slightly on the scaling of the signal. 
Therefore, it is best to normalize the NMR signal to a maximum of 1 
(this is default for the invert function, no need to do anything).


for convenience.

Note that for this method, the role of α is inverted, with larger α values
leading to less smoothing, not more. For that reason, `alpha=gcv()` (Mitchell method)
will not work. Please use `alpha=gcv(starting_value)` or `alpha=gcv(lower_limit, upper_limit)', or some lcurve method instead.
"""
struct pdhgm <: regularization_solver
    σ::Real
    τ::Real
    tol::Real
end
pdhgm(;sigma::Real=10, tau::Real=0.1, tol::Real=1e-5) = pdhgm(sigma, tau, tol)

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

    large = length(K) > 5000 ? true : false
    n = size(K, 2)

    Y, Ỹ, temp_vec = zeros.((n, n, n))
    f, f̃, f_prev  = ones.((n, n, n))
    
    Ks = (τ * α) .* (K' * s)   

    if large
        C = cholesky(I + (τ * α) .* (K' * K)) 
    else
        B = inv(I + τ * α * K' * K)
    end

    ε = tol + 1.0

    while ε > tol

        @. Ỹ = Y + σ * f̃
        @. Y = Ỹ / max(1, abs(Ỹ))
        
        @. temp_vec = f̃ - (τ * Y) + Ks
        large ? ldiv!(f, C, temp_vec) : mul!(f, B, temp_vec) 
        @. f = max(0, f)
        
        @. temp_vec = f - f_prev 
        ε = norm(temp_vec) / norm(f_prev)
        
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
