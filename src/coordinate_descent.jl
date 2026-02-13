# coordinate descent L1 regularisation, based on:
# https://github.com/RoyiAvital/Projects/blob/c7de61f03b1797d96d959b8b6ad718f5a4f040f8/Optimization/LsL1SolversAnalysis/SolveLsL1Cd.m

using LinearAlgebra

function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::cdL1)
    return solve_ls_l1_cd(K, g, α, solver.iterations, solver.tol)
end

# has been split into two functions, to make iter and tol type-stable
function solve_ls_l1_cd(
    K::AbstractMatrix{T}, 
    g::AbstractVector{T}, 
    α::Real, 
    iter::Int, 
    tol::Real
) where T

    n = size(K, 2)
    
    # Pre-calculate squared column norms and their reciprocals
    col_norms_sq = vec(sum(abs2, K, dims=1))
    inv_col_norms = 1 ./ col_norms_sq
    
    f = zeros(T, n)
    r = K*f - g

    # mX = zeros(T, n, iter)

    for i in 2:iter
        δ_max = zero(T) 
        
        for j in 1:n
            vK = view(K, :, j)
            old_fj = f[j]
            
            β = col_norms_sq[j] * old_fj - dot(vK, r)
            
            # Soft-thresholding
            new_fj = sign(β) * max(abs(β) - α, zero(T)) * inv_col_norms[j]

            δ = abs(new_fj - old_fj)
            if δ > δ_max
                δ_max = δ
            end           

            # Update residual: r_new = r_old + Kj * (new_fj - old_fj)
            diff = new_fj - old_fj
            if abs(diff) > 1e-15 # Optimization: only update if change is meaningful
                axpy!(diff, vK, r) 
            end
            
            f[j] = new_fj
        end

        # mX[:, i] .= f

        if δ_max < tol
            return f, r
        end
    end

    return f, r
end
