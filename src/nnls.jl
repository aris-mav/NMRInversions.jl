using LinearAlgebra, NonNegLeastSquares

export nnls
"""
    nnls(; L, algorithm)
Non-negative least squares method for regularization, 
implemented using the `NonNegLeastSquares` package extension.
Decent speed for small, 1D regularizations, but very slow for larger, 2D problems.
It can be used as a "solver" for invert function.

- `L` determines the tikhonov matrix.\
Can be either a matrix or an integer `n`. The latter will create a finite difference matrix\
of order `n`, which means that the penalty term will be the `n`'th derivative of the distribution.\
Defaults to `0`, where `L` becomes the Identity matrix.
- `algorithm` is one of:
    * `:nnls`
    * `:fnnls`
    * `:pivot`
"""
struct nnls <: regularization_solver 
    L::Union{Int, AbstractMatrix}
    algorithm::Symbol
end

nnls(;L=0, algorithm=:nnls) = nnls(L, algorithm)


function NMRInversions.solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::nnls)

    L = if isa(solver.L,Int)
        NMRInversions.Γ(size(K, 2), solver.L)
    else
        if size(solver.L, 2) == size(K, 2) 
            solver.L 
        else
            throw("Size mismatch between `solver.L` and the Kernel.")
        end
    end

    A = [K; √(α) .* L]
    b = [g; zeros(size(A, 1) - size(g, 1))]

    f = vec(nonneg_lsq(A, b ; alg=solver.algorithm))

    r = K * f - g

    return f, r
end
