"""
    sigmoid
using Base: return_types
Return an ascending sigmoid function of the same size as X.

Arguments:
- `X` is the x-axis of the data.
- `X_c` the value at which the sigmoid will be 0.

Keyword (optional) arguments:
- `width` determines the steepness of the sigmoid.
- `scale` function describing the scale of the input 
vector `X`. Defaults to `log10`, and `identity` can be 
used for a linear scale.

Some indicative values:
- `width=1`: almost like a line with slope `1`.
- `width=0.05`: looks like a typical sigmoid.
- `width=0.001`: makes it basically a step function.
"""
function sigmoid(X::AbstractVector, X_c::Real; 
                 width::Real=0.05, scale::Function=log10)

    scaled_X = scale.(X)
    sigcent = scaled_X .- scale(X_c)
    
    range_X = maximum(scaled_X) - minimum(scaled_X)
    effective_width = width * range_X
    
    return 1 ./ (1 .+ exp.(-sigcent ./ effective_width))

end


"""
    sge(;kwargs...)
Simultaneous Gaussian-Exponential inversion solver.
Wraps nnls() and utilises sigmoid functions as the diagonal 
of the Tikhonov matrix to penaltise gaussian components above the 
centre value, and exponential components below the centre value.

The following optional (Keyword) arguments can be provided.

- `s` is the sigmoid vector. If not provided, it will be constructed from the values below:
- `s_centre` default value is 1e-4 (seconds).
- `s_width` default value is 0.05.
- `s_weight` default value is 1.

"""
struct sge <: regularization_solver
    s::Vector{<:Real}
    s_centre::Real
    s_width::Real
    s_weight::Real
    algorithm::Symbol
end

sge(;s_centre::Real=1e-4, 
    s_width::Real=0.05, 
    s_weight::Real=0.001, 
    s::Vector{<:Real}=zeros(0),
    algorithm::Symbol = :nnls
    ) = sge(s, s_centre, s_width, s_weight, algorithm)

function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::sge)

    solver = if !isempty(solver.s)
        l1 = length(solver.s) * 2
        l2 = size(K, 2)
        if l1 == l2
            nnls(L = Diagonal(√α .+ [solver.s ; maximum(solver.s) .- solver.s]), algorithm = :nnls)
        else
            error("Length of sigmoid should be $l2 to match kernel size.")
        end
    else
        error("Provide a sigmoid in the sge solver.")
    end

    # using α=1 since we incorporate it in the sigmoid above
    f, r = solve_regularization(K, g, 1, solver)

end

function log_gaussian(x, μ, σ, amp)
    return amp .* exp.(-(log10.(x) .- log10(μ)).^2 ./ (2 * σ^2))
end

function test_sge()

    X = collect(logrange(1e-5, 1, 100))
    f_gauss = log_gaussian(X, 0.2e-3, 0.2, 0.6)
    f_exp = log_gaussian(X, 0.6e-1, 0.2, 0.4)

    x = collect(range(1e-5,1e-2,500))

    K_g = create_kernel(CPMG, x, X, gaussian = true)
    K_e = create_kernel(CPMG, x, X)


    y_gaussian = K_g * f_gauss
    y_exp = K_e * f_exp

    y = y_gaussian + y_exp

    data = input1D(CPMG, x, y)

    f = invert(data, solver=sge(), alpha=lcurve(1e-6,1e1), lims=X).f

end
