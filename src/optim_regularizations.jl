function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::brd)

    f = x -> objective_f(x, (α, g, K))
    g! = (grad, x) -> gradient_f!(grad, x, (α, g, K))
    h! = (hess, x) -> hessian_f!(hess, x, (α, g, K))

    c = Optim.optimize(
        f, g!, h!,
        ones(size(K, 1)),
        NewtonTrustRegion(),
        Optim.Options(x_tol = 1e-8)
    ).minimizer

    f = vec(max.(0, (K' * c)))
    r = K * f - g

    return f, r
end

function objective_f(c::AbstractVector, p::Tuple{<:Real, AbstractVector, AbstractMatrix } )
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0

    0.5 * c' * ((p[1] * I + (G * p[3]')) * c) - c' * p[2]

end

function gradient_f!(gradient, c, p::Tuple{<:Real, AbstractVector, AbstractMatrix })
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    gradient .= (p[1] * I + G * p[3]') * c - p[2]
end

function hessian_f!(H, c, p::Tuple{<:Real, AbstractVector, AbstractMatrix })
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    H .= p[1] * I + G * p[3]'
end


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
