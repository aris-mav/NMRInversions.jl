
export mexp
"""
    mexp(seq, u, x)
Construct a multi-exponential function from parameters u, and evaluate it at the values x. \n
The `u` vector should be of the form [a1, b1, a2, b2, ...], 
where a's are the amplitudes and b's are the reciprocals of the decay constants.
The length of `u` must be an even number.

Examples: \n
mexp(CPMG, [a,b] , x) = a * exp.( (1/b) * x) \n
mexp(CPMG, [a,b,c,d] , x) = a * exp.( (1/b) * x) + c * exp.((1/d) * x) \n
mexp(CPMG, [a,b,c,d,e,f] , x) = a * exp.( (1/b) * x) + c * exp.((1/d) * x) + e * exp.((1/f) * x) \n
(where a,b,c,d,e,f are all numbers)
. . .

"""
function mexp(seq, u, x)

    if seq == IR
        un = sum([u[i] for i in 1:2:length(u)])
        return un .- 2 .* mexp(u,x)

    elseif seq == SR
        un = sum([u[i] for i in 1:2:length(u)])
        return un .-  mexp(u,x)

    elseif seq == CPMG
        return mexp(u, x)

    elseif seq == PFG
        return mexpd(u, x)
    end

end
mexp(u, x) = sum([u[i] * exp.(-(1 / u[i+1]) * x) for i in 1:2:length(u)])
mexpd(u, x) = sum([u[i] * exp.(-(1 * u[i+1]) * x) for i in 1:2:length(u)])


"Loss function for multi-exponential fit, returns the sum of squared differences between the data and the model."
function mexp_loss(u, p)
    x = p[1]
    y = p[2]
    seq = p[3]
    return sum((mexp(seq, u, x) .- y) .^ 2)
end


export expfit
"""
    expfit(n, seq, x, y;solver)

Fit an exponential function to the data `x`, `y` using `n` exponential terms. \n
Starting points for the nonlinear regression are automatically chosen.
The outut is an `expfit_struct` structure.

Arguments:

- `n` : number of exponential terms.
- `seq` : pulse sequence.
- `x` : acquisition x parameter (time for relaxation or b-factor for diffusion).
- `y` : acquisition y parameter (magnetization).
- `solver` : optimization solver (default is BFGS).


Note that initial conditions are automatically determined. 
If you get NaN for any of the resulting parameters, 
try changing the initial conditions by calling the funtion using
the `u0` argument instead of `n`.
"""
function expfit(n::Int, seq::Type{<:NMRInversions.pulse_sequence1D}, x::Vector, y::Vector; kwargs...)

    # Make an educated guess of u0
    if seq == PFG
        u0 = [v for l in 1:n for v in ( abs(y[1])/(3*l -2) ,maximum(x)/100 * 10.0^(-(l/2 -0.5)) ) ]
    else 
        u0 = [v for l in 1:n for v in ( abs(y[1])/(3*l -2) ,maximum(x)/5 * 10.0^(-(l/2)) ) ]
    end

    return expfit(u0, seq, x, y; kwargs...)
end

"""
    expfit(n, data::input1D; kwargs...)
Similar to the `invert` fucntion, `expfit` can be called using an `input1D` structure.
"""
function expfit(u0::Vector{<:Real}, res::NMRInversions.input1D; kwargs...)
    expfit(u0, res.seq, res.x, res.y, kwargs...)
end


"""
    expfit(u0, seq, x, y; solver=BFGS())

Fit an exponential function to the data `x`, `y`. \n
The outut is an `expfit_struct` structure.

Arguments:

- `u0` : initial parameter guesses.
- `seq` : pulse sequence.
- `x` : acquisition x parameter (time for relaxation or b-factor for diffusion).
- `y` : acquisition y parameter (magnetization).
- `solver` : OptimizationOptimJL solver, defeault choice is BFGS().

The `u0` argument is a vector of initial parameter guesses, 
and it also determines the number of exponential terms used.
It should be of the form [a1, b1, a2, b2, ...], 
where a's are the amplitudes and b's are the reciprocals of the decay constants.
The length of `u0` must be an even number.

The following examples might help to clarify: \n
`expfit([a,b] , CPMG, x, y)` -> mono-exponential fit with initial guess: a * exp.( (1/b) * x) \n

`expfit([a,b,c,d] , CPMG, x, y)` -> double-exponential fit with initial guess: a * exp.( (1/b) * x) + c * exp.((1/d) * x) \n

`expfit([a,b,c,d,e,f] , CPMG, x, y)` -> triple-exponential fit with initial guess: a * exp.( (1/b) * x) + c * exp.((1/d) * x) + e * exp.((1/f) * x) \n

(where a,b,c,d,e,f are numbers of your choice)
Numbers of parameters beyond tri-exponential can also be used, but it is not recommended.

"""
function expfit(u0::Vector{<:Real}, seq::Type{<:NMRInversions.pulse_sequence1D}, x::Vector, y::Vector;
     solver=OptimizationOptimJL.BFGS())

    u0 = Float64.(u0)

    # Solve the optimization
    optf = Optimization.OptimizationFunction(mexp_loss, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(optf, u0, (x, y, seq), lb=zeros(length(u0)), ub=Inf .* ones(length(u0)))
    u = OptimizationOptimJL.solve(prob, solver, maxiters=5000, maxtime=100)

    # Determine what's the x-axis of the seq (time or bfactor)
    seq == NMRInversions.PFG ? x_ax = "b" : x_ax = "t"

    # Determine whether there should be a multiplication or a division in the exp() formula 
    seq == NMRInversions.PFG ?  op = "*" : op = "/"

    un = sum([u[i] for i in 1:2:length(u)]) 

    # Get the fit's equation as a string
    eq = join( [(i == 1 ? "" : " + ") * string(round(u[i], sigdigits=2)) * " * exp(-$x_ax" * op * string(round(u[i+1], sigdigits=2)) * ")" for i in 1:2:length(u)])

    # Normalised version of the fit equation
    eqn = string(round(un, sigdigits=2)) * " * (" * join([(i == 1 ? "" : " + ") * string(round(u[i] / un, sigdigits=2)) * " * exp(-$x_ax" * op * string(round(u[i+1], sigdigits=2)) * ")" for i in 1:2:length(u)]) * ")"

    if seq == IR
        unstr = string(round(un,sigdigits = 2))
        eq = unstr * " - 2 * (" * eq * ")" 
        eqn = unstr * " - 2 * " *  eqn   

    elseif seq == SR
        unstr = string(round(un,sigdigits = 2))
        eq = unstr * " -" * "(" * eq * ")" 
        eqn = unstr * " -" *  eqn   
    end

    # Calculate the residuals
    r = mexp(seq, u, x) - y

    if length(u0) == 2
        return expfit_struct(seq, x, y, u, u0, r, eq, eq, "")
    else
        return expfit_struct(seq, x, y, u, u0, r, eq, eqn, "")
    end

end

"""
    expfit(u0, data::input1D; kwargs...)
Similar to the `invert` fucntion, `expfit` can be called using an `input1D` structure.
"""
function expfit(n::Int, res::NMRInversions.input1D; kwargs...)
    expfit(n, res.seq, res.x, res.y, kwargs...)
end

