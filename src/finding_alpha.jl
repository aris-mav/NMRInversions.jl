export gcv, lcurve


## Mitchell et al method
"""
    gcv_mitchell
Generalized cross validation for finding the optimal regularization parameter α,
following the method in Mitchell 2012.
"""
struct gcv_mitchell <: alpha_optimizer end
"""
    gcv()
Constructor for the gcv method described in Mitchell 2012.
No arguments required.
"""
gcv() = gcv_mitchell()



# UNIVARIATE METHODS
"""
Finding the optimal regularization parameter α,
using univariate optimization algorithm. Brent() is the default option, 
GoldenSection() can be used as an alternative.
"""
struct find_alpha_univariate <: alpha_optimizer
    search_method::Symbol
    lower::Real
    upper::Real
    algorithm::Optim.UnivariateOptimizer
    abs_tol::Real
    rel_tol::Real
end

"""
    gcv(lower, upper ; kwargs...)
Constructor for finding the optimal alpha value via gcv score
univariate optimization, given some lower and upper bounds.

 Necessary (positional) arguments:
- `lower` is the lower bound, or lowest alpha value the algorighm will consider.
- `upper` is the upper bound, or highest alpha value the algorighm will consider.

 Optional (keyword) arguments:
- `algorithm` determines which method will be used by Optim.jl to solve the problem.
Default is `Brent()` and the only alternative is `GoldenSection()`, which is normally slower.
- `abs_tol` determines what's the smallest absolute change in the gcv score the algorithm should care about.
Default is `1e-3`.
- `rel_tol` determines what's the smallest relative change in the gcv score the algorithm should care about.
Default is `1e-1`.
"""
gcv(lower::Real, upper::Real; algorithm = Brent(), abs_tol=1e-3, rel_tol=1e-1) = 
    find_alpha_univariate(:gcv, Float64(lower), Float64(upper), algorithm, abs_tol, rel_tol)

"""
    lcurve(lower, upper ; kwargs...)
Constructor for finding the optimal alpha value via lcurve curvature
univariate optimization, given some lower and upper bounds.

 Necessary (positional) arguments:
- `lower` is the lower bound, or lowest alpha value the algorighm will consider.
- `upper` is the upper bound, or highest alpha value the algorighm will consider.

 Optional (keyword) arguments:
- `algorithm` determines which method will be used by Optim.jl to solve the problem.
Default is `Brent()` and the only alternative is `GoldenSection()`, which is normally slower.
- `abs_tol` determines what's the smallest absolute change in the lcurve score the algorithm should care about.
Default is `1e-3`.
- `rel_tol` determines what's the smallest relative change in the lcurve score the algorithm should care about.
Default is `1e-1`.
"""
lcurve(lower::Real, upper::Real; algorithm = Brent(), abs_tol=1e-3, rel_tol=1e-1) = 
    find_alpha_univariate(:lcurve, Float64(lower), Float64(upper), algorithm, abs_tol, rel_tol)


## Fminbox methods

"""
Finding the optimal regularization parameter α,
using Fminbox optimization algorithm (LBFGS by
default). 
"""
struct find_alpha_box <: alpha_optimizer
    search_method::Symbol
    start::Real
    algorithm::Optim.FirstOrderOptimizer
    opts::Optim.Options
end
"""
    gcv(start; kwargs...)
Constructor for finding the optimal alpha value via gcv score
box optimization, given a starting value.

 Necessary (positional) arguments:
- `start` is the starting alpha value. 
Choose something sensible, usually a value between 0.1 and 10 would work well.

 Optional (keyword) arguments:
- `algorithm` determines which method will be used by Optim.jl to solve the problem.
Default is LBFGS(). Only first order optimizers can be chosen here. 
For more details, refer to Optim.jl documentation.
- `opts` an Optim.Options() structure which can provide some preferences to the solver.
"""
gcv(start::Real; algorithm = LBFGS(), opts = Optim.Options(x_abstol=1e-3)) = 
    find_alpha_box(:gcv, Float64(start), algorithm, opts)

"""
    lcurve(start; kwargs...)
Constructor for finding the optimal alpha value via lcurve curvature
box optimization, given a starting value.

 Necessary (positional) arguments:
- `start` is the starting alpha value. 
Choose something sensible, usually a value between 0.1 and 10 would work well.

 Optional (keyword) arguments:
- `algorithm` determines which method will be used by Optim.jl to solve the problem.
Default is LBFGS(). Only first order optimizers can be chosen here. 
For more details, refer to Optim.jl documentation.
- `opts` an Optim.Options() structure which can provide some preferences to the solver.
Please have a look [here](https://julianlsolvers.github.io/Optim.jl/v1.10/user/config/).
"""
lcurve(start::Real; algorithm = LBFGS(), opts = Optim.Options(x_abstol = 1e-3)) = 
    find_alpha_box(:lcurve, Float64(start), algorithm, opts)




# Full lcurve method

struct lcurve_range <: alpha_optimizer 
    lowest_value::Real
    highest_value::Real
    number_of_steps::Int
end
"""
    lcurve(lowest_value, highest_value, number_of_steps)
Constructor for testing all L-curve curvatures between some bounds and 
picking the optimal.

- `lowest_value` is the lowest alpha value.
- `highest_value` is the highest alpha value.
- `number_of_steps` is the number of alpha values that will be tested.

This is a very crude and rather slow method, mostly for demonstration purposes.
"""
lcurve(lowest_value::Real, highest_value::Real, number_of_steps::Int,) = 
    lcurve_range(lowest_value::Real, highest_value::Real, number_of_steps::Int)



# Functions


"""
Compute the Generalized Cross Validation (GCV) score for a regularization solution, and return the score, as well as the next α to test. 
(look Mitchell 2012 paper https://doi.org/10.1016/j.pnmrs.2011.07.002)

- `α` : smoothing term
- `r` : residuals from the regularization (`r = Kf .- g`)
- `s` : singular values of the kernel
- `x` : `Ṽ₀'f` or `s .* c`

The symbols above follow the notation used in the paper.

"""
function gcv_score(α, r, s, x; next_alpha=true)

    ñ = length(s)
    c = r ./ α
    σ² = s .^ 2
    m̂ = sum(σ² ./ (σ² .+ α))
    φₙ = (α .^ 2 * c' * c .* ñ) ./ ((ñ - m̂) .^ 2)  # GCV score to be minimized

    if next_alpha

        dm̂ = sum(σ² ./ ((σ² .+ α) .^ 2))
        t̂ = sum(x .^ 2 ./ (σ² .+ α))
        αₙ = (α .^ 2 * c' * c * dm̂) / (t̂ * (ñ - m̂))  # α update value, test this one next
        return φₙ, αₙ

    else
        return φₙ
    end
end

function gcv_cost(α::Real, 
                  svds::svd_kernel_struct, 
                  solver::regularization_solver)

    f, r = solve_regularization(svds.K, svds.g, α, solver)
    return gcv_score(α, r, svds.S, (svds.V' * f), next_alpha = false) 
end

"""
Compute the curvature of the L-curve at a given point.
(Hansen 2010 page 92-93)

- `f` : solution vector
- `r` : residuals
- `α` : smoothing term
- `A` : Augmented kernel matrix (`K` and `αI` stacked vertically)
- `b` : Augmented residuals (`r` and `0` stacked vertically)

"""
function l_curvature(f, r, α, A, b)

    ξ = dot(f, f)
    ρ = dot(r, r)
    λ = sqrt(α)

    z = A \ b
    ∂ξ∂λ = (4 / λ) * dot(f, z)

    ĉ = 2 * (ξ * ρ / ∂ξ∂λ) * (α * ∂ξ∂λ * ρ + 2 * ξ * λ * ρ + λ^4 * ξ * ∂ξ∂λ) / ((α * ξ^2 + ρ^2)^(3 / 2))

    return ĉ
end

function l_cost(K, g, α, solver)

    f, r = NMRInversions.solve_regularization(K, g, α, solver)

    A = sparse([K; √(α) * LinearAlgebra.I ])
    b = [r; zeros(size(A, 1) - size(r, 1))]

    return l_curvature(f, r, α, A, b)

end


"""
Solve repeatedly until the GCV score stops decreasing, following Mitchell 2012 paper.
Select the solution with minimum gcv score and return it, along with the residuals.
"""
function find_alpha(svds::svd_kernel_struct, 
                   solver::regularization_solver,
                   mode::gcv_mitchell ; silent = false)

    s̃ = svds.S
    ñ = length(s̃)

    #Initial guess (overestimate)
    α = sum(s̃ .^ 2) / ñ

    if any(isa.( (solver,), (
        cdL1,
        pdhgm,
    )))
        throw("The Mitchell et. al. gcv method will not work a `$(typeof(solver))` solver. \
              Please use a different alpha option, or change the solver \
              (e.g.: `alpha=gcv(0.01,10)` or `solver=brd()`).")
    end

    alphas = []
    phis = []
    f_star = []

    conditions = [true for _ in 1:3]

    while all(conditions) 

        silent ? nothing : display("Testing α = $(round(α,sigdigits=3))")

        f, r = solve_regularization(svds.K, svds.g, α, solver)
        push!(alphas, α) # Add the just tested α to the array

        φ, α = gcv_score(α, r, svds.S, (svds.V' * f)) # Compute φ for current α, and also compute new α 
        push!(phis, φ)

        conditions .= [
            length(phis) > 1 ? phis[end] < phis[end-1] : true,
            length(alphas) > 1 ? (abs(alphas[end] - alphas[end-1]) > 1e-3) : true,
            length(alphas) < 15
        ]

        # If the gcv score dropped, save the solution
        if conditions[1] 
            f_star = copy(f)
        end

    end

    silent ? nothing : display("The optimal α is $(round(alphas[end-1], sigdigits=3))")
    α = alphas[end-1]
    f = f_star
    r = svds.g - svds.K * f

    return f, r, α
end


"""
Find alpha via univariate optimization.
"""
function find_alpha(svds::svd_kernel_struct,
                   solver::regularization_solver,
                   mode::find_alpha_univariate; silent = true
                   )

    local f
    if mode.search_method == :gcv
        f = x -> gcv_cost(x, svds, solver)

    elseif mode.search_method == :lcurve
        f = x -> l_cost(svds.K, svds.g, x, solver)

    end

    sol = optimize(
        f,
        mode.lower, mode.upper,
        mode.algorithm,
        abs_tol = mode.abs_tol,
        rel_tol = mode.rel_tol,
        show_trace = !silent
    )

    α = sol.minimizer
    silent ? nothing : display("Converged at α =$(round(α,sigdigits=3)), after $(sol.f_calls) calls.")

    f, r = NMRInversions.solve_regularization(svds.K, svds.g, α, solver)

    return f, r, α

end

 
"""
Find alpha using Fminbox optimization.
"""
function find_alpha(svds::svd_kernel_struct,
                    solver::regularization_solver,
                    mode::find_alpha_box ; silent = true
                    )

    local f
    if mode.search_method == :gcv
        f = x -> gcv_cost(first(x), svds, solver)

    elseif mode.search_method == :lcurve

        f = x -> l_cost(svds.K, svds.g, first(x), solver)
    end

    sol = optimize(
        f,
        [0], [Inf], [mode.start],
        Fminbox(mode.algorithm),
        change_field(mode.opts, show_trace = !silent)
    )

    α = sol.minimizer[1]

    silent ? nothing : display("Converged at α =$(round(α,sigdigits=3)), after $(sol.f_calls) calls.")
    f, r = NMRInversions.solve_regularization(svds.K, svds.g, α, solver)

    return f, r, α

end


"""
Test `n` alpha values between `lower` and `upper` and select the one 
which is at the heel of the L curve, accoding to Hansen 2010.
"""
function find_alpha(svds::svd_kernel_struct,
                   solver::regularization_solver,
                   mode::lcurve_range; silent = false
                   )

    alphas = exp10.(range(log10(mode.lowest_value), 
                          log10(mode.highest_value), 
                          mode.number_of_steps))

    curvatures = zeros(length(alphas))

    for (i, α) in enumerate(alphas)

        curvatures[i] = l_cost(svds.K, svds.g, α, solver)

    end

    α = alphas[argmin(curvatures)]
    silent ? nothing : display("The optimal α is $(round(α,sigdigits=3))")

    f, r = NMRInversions.solve_regularization(svds.K, svds.g, α, solver)
    return f, r, α

end
