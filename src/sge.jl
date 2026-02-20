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


# Testing the paper method, for reference
function SGE_washburn(X, M, to, T2;
                      smoothing = 0.0001,
                      sig = 38,
                      sigwidth = 1.0,
                      sigweightA = 0.001,
                      sigweightB = 0.001,
                      )   

    halfX = div(length(X), 2)
    A = X[1:halfX]
    B = X[(halfX + 1):end]
    
    T = eltype(X) 
    
    Mhat = zeros(T, length(M))
    
    sigcent = collect(T, (-sig):(-sig + halfX - 1))
    
    sigup = ((1 ./ (1 .+ exp.(-sigcent .* sigwidth))) .* sigweightA) .+ smoothing
    sigdown = (smoothing .+ (1 .- (1 ./ (1 .+ exp.(-sigcent .* sigwidth)))) .* sigweightB)
    
    for ii in 1:halfX
        Mhat .+= A[ii] .* exp.(-((to ./ T2[ii]) .^ 2)) .+ B[ii] .* exp.(-to ./ T2[ii])
    end
    
    SS = sum((M .- Mhat) .^ 2) + sum(A .* sigup) + sum(B .* sigdown)
    
    return SS
end

function log_gaussian(x, μ, σ, amp)
    return amp .* exp.(-(log10.(x) .- log10(μ)).^2 ./ (2 * σ^2))
end

function test_sge()

    # make T distributions
    X = collect(logrange(1e-5, 1, 50))
    f_gauss = log_gaussian(X, 0.2e-3, 0.2, 0.6)
    f_exp = log_gaussian(X, 0.6e-1, 0.2, 0.4)

    x = collect(range(1e-5,1e-2,500))

    K_g = create_kernel(CPMG, x, X, gaussian = true)
    K_e = create_kernel(CPMG, x, X)

    y_gaussian = K_g * f_gauss
    y_exp = K_e * f_exp

    y = y_gaussian + y_exp

    s = sigmoid(X, 5e-3 , width = 0.001)

    # return X, f_gauss + f_exp

    K_sge = [K_g K_e]

    s_weight = 1

    solver = optim_nnls(L = Diagonal(s_weight .* [s ; 1 .- s]), opts= Optim.Options(show_trace=true))

    f, r = solve_regularization(K_sge, y, 1.0, solver)

    return f

    # initial_guess = fill(1.0 / (length(X) * 2), length(X) * 2)
    # lower_bounds = fill(0.0, length(initial_guess))
    # upper_bounds = fill(200.0, length(initial_guess))
    #
    # inner_obj = X -> NMRInversions.SGE_washburn(X, y, x, X)
    #
    # res = optimize(
    #     inner_obj, 
    #     lower_bounds, 
    #     upper_bounds, 
    #     initial_guess, 
    #     Fminbox(LBFGS()), 
    #     Optim.Options(
    #         show_trace = true,
    #         outer_iterations = 5,
    #         iterations = 100,
    #         f_abstol = 1e-3,
    #         f_reltol = 1e-3,
    #         g_tol = 1e-4
    #     ),
    #     autodiff = :forward
    # )
    #
    # optimized_X = Optim.minimizer(res)
    # A_final = optimized_X[1:length(X)]
    # B_final = optimized_X[length(X)+1:end]
    #
    # return x, y, X, f_gauss, f_exp, A_final, B_final
    #
end

#==
function testplot(x, y, T2Values, f_gauss, f_exp, A_final, B_final)
    fig = GLMakie.Figure()
    ax_1 = GLMakie.Axis(fig[1,:], xscale = log10, title="Distribution")
    ax_2 = GLMakie.Axis(fig[2,:], title = "CPMG")
    ax_3 = GLMakie.Axis(fig[3,:], xscale = log10, title="Results")
    lines!(ax_1, T2Values, f_gauss)
    lines!(ax_1, T2Values, f_exp)
    lines!(ax_2, x, y)
    lines!(ax_3, T2Values, A_final)
    lines!(ax_3, T2Values, B_final)
    return fig
end
==#
