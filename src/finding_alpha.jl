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
                  solver::Union{regularization_solver, Type{<:regularization_solver}})

    #=display("Testing α = $(round(α,sigdigits=3))")=#
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

    ξ = f'f
    ρ = r'r
    λ = √α

    z = NMRInversions.solve_ls(A, b)

    ∂ξ∂λ = (4 / λ) * f'z

    ĉ = 2 * (ξ * ρ / ∂ξ∂λ) * (α * ∂ξ∂λ * ρ + 2 * ξ * λ * ρ + λ^4 * ξ * ∂ξ∂λ) / ((α * ξ^2 + ρ^2)^(3 / 2))

    return ĉ

end

function l_cost(K, g, α, solver)

    f, r = NMRInversions.solve_regularization(K, g, α, solver)

    A = sparse([K; √(α) * LinearAlgebra.I ])
    b = sparse([r; zeros(size(A, 1) - size(r, 1))])

    return l_curvature(f, r, α, A, b)

end


"""
Solve repeatedly until the GCV score stops decreasing, following Mitchell 2012 paper.
Select the solution with minimum gcv score and return it, along with the residuals.
"""
function find_alpha(svds::svd_kernel_struct, 
                   solver::Union{regularization_solver, Type{<:regularization_solver}},
                   mode::gcv_mitchell ; silent = false)

    s̃ = svds.S
    ñ = length(s̃)

    #Initial guess (overestimate)
    α = sum(s̃ .^ 2) / ñ

    if any(isa.( (solver,), (
        cdL1,
        pdhgm,
    )))
        throw("The Mitchell et. al. gcv method will not work a `$(typeof(solver)) solver. \
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
                   solver::Union{regularization_solver, Type{<:regularization_solver}},
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
        show_trace = !silent
    )

    α = sol.minimizer
    silent ? nothing : display("Converged at α =$(round(α,sigdigits=3)), after $(sol.f_calls) calls.")

    f, r = NMRInversions.solve_regularization(svds.K, svds.g, α, solver)

    return f, r, α

end

 
# used below to add silend option to Optim.Options immutable struct
function change_field(x::T; kwargs...) where {T}
    # Convert keyword overrides to a dictionary for fast lookup
    kw = Dict(kwargs)

    # Get the fields in the correct positional order
    vals = [ get(kw, f, getfield(x, f)) for f in fieldnames(T) ]

    # Call the canonical positional constructor
    return T(vals...)
end


"""
Find alpha using Fminbox optimization.
"""
function find_alpha(svds::svd_kernel_struct,
                    solver::Union{regularization_solver, Type{<:regularization_solver}},
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
                   solver::Union{regularization_solver, Type{<:regularization_solver}},
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
