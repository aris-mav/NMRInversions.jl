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


"""
Solve repeatedly until the GCV score stops decreasing.
Select the solution with minimum gcv score and return it, along with the residuals.
"""
function solve_gcv(svds::svd_kernel_struct, solver::Union{regularization_solver, Type{<:regularization_solver}})

    s̃ = svds.S
    ñ = length(s̃)

    #Initial guess (overestimate)
    α = sum(s̃ .^ 2) / ñ

    alphas = []
    phis = []
    f_star = []

    conditions = [true for _ in 1:3]

    while all(conditions) 

        display("Testing α = $(round(α,sigdigits=3))")

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

    display("The optimal α is $(round(alphas[end-1], digits=3))")
    α = alphas[end-1]
    f = f_star
    r = svds.g - svds.K * f

    return f, r, α
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


"""
Test `n` alpha values between `lower` and `upper` and select the one 
which is at the heel of the L curve, accoding to Hansen 2010.
"""
function solve_l_curve(K, g, solver, lower, upper, n)

    alphas = exp10.(range(log10(lower), log10(upper), n))
    curvatures = zeros(length(alphas))

    for (i, α) in enumerate(alphas)
        display("Testing α = $(round(α,sigdigits=3))")

        f, r = NMRInversions.solve_regularization(K, g, α, solver)

        A = sparse([K; √(α) * LinearAlgebra.I ])
        b = sparse([r; zeros(size(A, 1) - size(r, 1))])

        c = l_curvature(f, r, α, A, b)

        curvatures[i] = c

    end

    α = alphas[argmin(curvatures)]
    display("The optimal α is $(round(α,sigdigits=3))")

    f, r = NMRInversions.solve_regularization(K, g, α, solver)

    return f, r, α

end
