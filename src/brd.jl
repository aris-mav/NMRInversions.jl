function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::brd)

    f = x -> objective_f(x, (α, g, K))
    g! = (grad, x) -> gradient_f!(grad, x, (α, g, K))
    h! = (hess, x) -> hessian_f!(hess, x, (α, g, K))

    c = Optim.optimize(
        f, g!, h!,
        ones(size(K, 1)),
        NewtonTrustRegion(),
        Optim.Options(x_abstol = 1e-8)
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
