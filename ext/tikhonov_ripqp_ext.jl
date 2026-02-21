module tikhonov_ripqp_ext

using NMRInversions
using RipQP
using QuadraticModels
using SparseArrays
using LinearAlgebra

struct ripqp <: regularization_solver end

function NMRInversions.solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::Type{NMRInversions.ripqp})

    ## solve 0.5 xᵀ H x + cᵀx + c₀ 
    ## s.t. Ax = b
    ## r = Kf-g
    ## x = [f;r]
    ## [K I]x = [g] => Kf = g + r

    b = g
    m, n = size(K)
    A = [SparseArrays.sparse(K) LinearAlgebra.I]
    c = zeros(m + n)
    lvar = [zeros(n); fill(-Inf, m)] # No bounds to residuals, f positive

    H = [2*α*I SparseArrays.spzeros(n, m);
        SparseArrays.spzeros(m, n) 2*LinearAlgebra.I]

    qm = QuadraticModel(c, H;
        A=A,
        lcon=b, ucon=b,
        lvar=lvar)

    stats = RipQP.ripqp(qm, display=false)

    return stats.solution[1:n], stats.solution[n+1:n+m]  # f, r

end




end
