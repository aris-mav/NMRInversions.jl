function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::optim_nnls)

    A = sparse([K; √(α) .* NMRInversions.Γ(size(K, 2), solver.order)])
    b = sparse([g; zeros(size(A, 1) - size(g, 1))])

    f = solve_nnls(A, b, L = solver.L)

    r = K * f - g

    return f, r
end


"""
Solve a least squares problem, with nonnegativity constraints.
"""
function solve_nnls(A::AbstractMatrix, b::AbstractVector ;
                    start = ones(size(A, 2)),
                    L = 2)

    x = Optim.optimize(
        x -> norm( A*x - b, L), 
        zeros(size(A, 2)), fill(Inf, size(A, 2)), start,
        autodiff = :forward,
    ).minimizer

    return x
end


"""
Solve a least squares problem, without nonnegativity constraints.
"""
function solve_ls(A::AbstractMatrix, b::AbstractVector)

    x = optimize(
        x -> norm( A*x - b, 2), 
        fill(-Inf, size(A, 2)), fill(Inf, size(A, 2)),
        ones(size(A, 2)),
        autodiff = :forward
    ).minimizer

    return x
end
