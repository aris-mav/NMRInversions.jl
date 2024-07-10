
struct brd_solver <: tikhonov_solver end
brd = brd_solver()
export brd

function minfun(c::AbstractVector, p=(α, data, K))
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    f = 0.5 * c' * ((p[1] * I + (G * p[3]')) * c) - c' * p[2]
end

function grad(gradient, c, p=(α, data, K))
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    gradient .= (p[1] * I + G * p[3]') * c - p[2]
end

function hess(H, c, p=(α, data, K))
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    H .= p[1] * I + G * p[3]'
end


function solve_tikhonov(K::AbstractMatrix, g::AbstractVector, α::Real, solver::NMRInversions.brd_solver, order::Int=0)

    optf = OptimizationFunction(minfun, grad=grad, hess=hess)
    prob = OptimizationProblem(optf, ones(size(K, 1)), (α, g, K))
    c = OptimizationOptimJL.solve(prob, NewtonTrustRegion(), x_tol=1e-8, maxiters=5000, maxtime=100)
    f = vec(max.(0, (K' * c)))
    r = K * f - g

    return f, r

    # if solver == :ripqp
    #     ## solve 0.5 xᵀ H x + cᵀx + c₀ 
    #     ## s.t. Ax = b
    #     ## r = Kf-g
    #     ## x = [f;r]
    #     ## [K I]x = [g] => Kf = g + r

    #     b = g
    #     m, n = size(K)
    #     A = [sparse(K) I]
    #     c = zeros(m + n)
    #     lvar = [zeros(n); fill(-Inf, m)] # No bounds to residuals, f positive

    #     H = [2*α*I spzeros(n, m);
    #         spzeros(m, n) 2*I]

    #     qm = QuadraticModel(c, H;
    #         A=A,
    #         lcon=b, ucon=b)
    #     # lvar=lvar)

    #     stats = ripqp(qm, display=false)

    #     return stats.solution[1:n], stats.solution[n+1:n+m]  # f, r

    # elseif solver == :brd

    #     optf = OptimizationFunction(minfun, grad=grad, hess=hess)
    #     prob = OptimizationProblem(optf, ones(size(K, 1)), (α, g, K))
    #     c = OptimizationOptimJL.solve(prob, NewtonTrustRegion(), x_tol=1e-8, maxiters=5000, maxtime=100)
    #     f = vec(max.(0, (K' * c)))
    #     r = K * f - g

    #     return f, r

    # elseif solver == :HiGHSlp

    #     A = sparse([K; √(α) .* Γ(size(K, 2), order)])
    #     b = sparse([g; zeros(size(A, 1) - size(g, 1))])

    #     model = Model(HiGHS.Optimizer)
    #     set_silent(model)
    #     @variable(model, f[1:size(K)[2]] >= 0)
    #     @objective(model, Min, sum((A * f .- b) .^ 2))

    #     optimize!(model)
    #     @assert is_solved_and_feasible(model)

    #     return value.(f), vec(K * value.(f) .- g)



    # else # other qp formulations

    #     # if solver == :Gurobi
    #     #     model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    #     # elseif solver == :HiGHS
    #     #     model = Model(HiGHS.Optimizer)
    #     # else
    #     #     error("Unknown solver")
    #     # end

    #     # set_silent(model)
    #     # @variable(model, f[1:size(K, 2)] >= 0)
    #     # @variable(model, r[1:size(K, 1)])
    #     # @constraint(model, r .== g .- K * f)
    #     # @objective(model, Min, LinearAlgebra.dot(r, r) + α * LinearAlgebra.dot(f, f))
    #     # optimize!(model)
    #     # @assert is_solved_and_feasible(model)
    #     # return value.(f), value.(r)
    # end

end