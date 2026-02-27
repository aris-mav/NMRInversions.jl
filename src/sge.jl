using Optim

"""
    sigmoid
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
    sge
Simultaneous Gaussian-Exponential inversion solver.
"""
struct sge <: regularization_solver
    s_centre::Real
    s_width::Real
    s_weight::Real
end

function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::sge)

    s = sigmoid(X, solver.s_centre , width = solver.s_width)

    solver = optim_nnls(L = Diagonal([s ; 1 .- s]))
    f, r = solve_regularization(K_sge, g, α, solver)

end

function log_gaussian(x, μ, σ, amp)
    return amp .* exp.(-(log10.(x) .- log10(μ)).^2 ./ (2 * σ^2))
end

function test_sge(mode::Int = 1)

    X = collect(logrange(1e-5, 1, floor(Int, 128/mode)))
    f_gauss = log_gaussian(X, 0.2e-3, 0.2, 0.6)
    f_exp = log_gaussian(X, 0.6e-1, 0.2, 0.4)

    x = collect(range(1e-5,1e-2,500))

    K_g = create_kernel(CPMG, x, X, gaussian = true)
    K_e = create_kernel(CPMG, x, X)

    y_gaussian = K_g * f_gauss
    y_exp = K_e * f_exp

    y = y_gaussian + y_exp

    s = sigmoid(X, 5e-3 , width = 0.001)

    K_sge = [K_g K_e]

    s_weight = 1

    if mode == 1

        # solver=optim_nnls(opts= Optim.Options(show_trace=true))
        # return @time invert(CPMG,x,y, alpha = 1.0, solver=solver )

        solver = nnls(L = 0, algorithm=:fnnls)
        return @time invert(CPMG,x,y,  solver=solver )

    elseif mode ==2

        # solver = optim_nnls(L = Diagonal(s_weight .* [s ; 1 .- s]), opts= Optim.Options(show_trace=true))
        solver = nnls(L = Diagonal(s_weight .* [s ; 1 .- s]),algorithm=:pivot)

        f, r = @time solve_regularization(K_sge, y, 0.1, solver)
        return f
    end
end
