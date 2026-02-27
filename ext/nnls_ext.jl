module nnls_ext

using NMRInversions
using LinearAlgebra, NonNegLeastSquares

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

end
