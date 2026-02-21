using Optim

export optim_nnls
"""
    optim_nnls(order)
Simple non-negative least squares method for regularization, 
implemented using Optim.jl.
All around effective, but can be slow for large problems, such as 2D inversions.
It can be used as a "solver" for invert function.

- `L` determines the tikhonov matrix.\
Can be either a matrix or an integer `n`. The latter will create a finite difference matrix\
of order `n`, which means that the penalty term will be the `n`'th derivative of the distribution.\
Defaults to `0`, where `L` becomes the Identity matrix.
- `algorithm` is of `Optim.FirstOrderOptimizer` type (defaults to `LBFGS()`).
- `opts` an `Optim.Options()` structure which can provide some preferences to the solver.
"""
struct optim_nnls <: regularization_solver 
    L::Union{Int, AbstractMatrix}
    algorithm::Optim.FirstOrderOptimizer
    opts::Optim.Options
end

optim_nnls(;L= 0, algorithm= Optim.LBFGS(), opts= Optim.Options()) = optim_nnls(L, algorithm, opts)

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
