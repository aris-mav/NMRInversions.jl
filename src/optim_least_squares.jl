function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::optim_nnls)

    L = if isa(solver.L,Int)
        NMRInversions.Γ(size(K, 2), solver.L)
    else
        if size(solver.L, 2) == size(K, 2) 
            solver.L 
        else
            throw("Size mismatch between `optim_nnls.L` and the Kernel.")
        end
    end

    A = sparse([K; √(α) .* L])
    b = sparse([g; zeros(size(A, 1) - size(g, 1))])

    f = solve_nnls(A, b, solver)

    r = K * f - g

    return f, r
end


"""
Solve a least squares problem, with nonnegativity constraints.
"""
function solve_nnls(A::AbstractMatrix, b::AbstractVector, solver::optim_nnls ;
                    start = ones(size(A, 2)))

    x = Optim.optimize(
        x -> norm(A*x - b), 
        zeros(size(A, 2)), fill(Inf, size(A, 2)), start,
        Fminbox(solver.algorithm),
        autodiff = :forward,
        solver.opts
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
