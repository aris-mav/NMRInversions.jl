
function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::Type{brd})

    optf = Optimization.OptimizationFunction(objective_f, grad=gradient_f, hess=hessian_f)
    prob = Optimization.OptimizationProblem(optf, ones(size(K, 1)), (α, g, K))
    c = OptimizationOptimJL.solve(prob, OptimizationOptimJL.NewtonTrustRegion(), x_tol=1e-8, maxiters=5000, maxtime=100)
    f = vec(max.(0, (K' * c)))
    r = K * f - g

    return f, r

end

function objective_f(c::AbstractVector, p=(α, data, K))
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    f = 0.5 * c' * ((p[1] * I + (G * p[3]')) * c) - c' * p[2]
end
function gradient_f(gradient, c, p=(α, data, K))
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    gradient .= (p[1] * I + G * p[3]') * c - p[2]
end
function hessian_f(H, c, p=(α, data, K))
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    H .= p[1] * I + G * p[3]'
end


function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::optim_nnls)

    A = sparse([K; √(α) .* NMRInversions.Γ(size(K, 2), solver.order)])
    b = sparse([g; zeros(size(A, 1) - size(g, 1))])

    f = solve_nnls(A, b)

    r = K * f - g

    return f, r
end


"""
Solve a least squares problem, with nonnegativity constraints.
"""
function solve_nnls(A::AbstractMatrix, b::AbstractVector)

    optf = Optimization.OptimizationFunction(obj_f, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(optf, ones(size(A, 2)), (A, b), lb=zeros(size(A, 2)), ub=Inf * ones(size(A, 2)))
    x = OptimizationOptimJL.solve(prob, OptimizationOptimJL.LBFGS(), maxiters=5000, maxtime=100)

    return x
end
function obj_f(x, p)
    return sum((p[1] * x - p[2]).^2)
end


"""
Solve a least squares problem, without nonnegativity constraints.
"""
function solve_ls(A::AbstractMatrix, b::AbstractVector)

    optf = Optimization.OptimizationFunction(obj_ls, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(optf, ones(size(A, 2)), (A, b), 
                                            lb= -Inf .* ones(size(A, 2)), ub=Inf .* ones(size(A, 2)))
    x = OptimizationOptimJL.solve(prob, OptimizationOptimJL.LBFGS(), maxiters=5000, maxtime=100)

    return x
end
function obj_ls(x, p)
    return sum((p[1] * x - p[2]).^2)
end
