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
function mexp(u::Vector{<:Number}, ax::DataAxis{1})

    if length(u) % 2 != 0
        throw(
            ErrorException("The length of `u` must be a multiple of 2.")
        )
    end

    if ax isa IR
        un = sum([u[i] for i in 1:2:length(u)])
        return un .- 2 .* mexp(u, values(ax))

    elseif ax isa SR
        un = sum([u[i] for i in 1:2:length(u)])
        return un .-  mexp(u, values(ax))

    elseif ax isa CPMG
        return mexp(u, values(ax))

    elseif ax isa PFG
        return mexpd(u, values(ax))
    end

end
mexp(u, x) = sum([u[i] * exp.(-(1 / u[i+1]) * x) for i in 1:2:length(u)])
mexpd(u, x) = sum([u[i] * exp.(-(1 * u[i+1]) * x) for i in 1:2:length(u)])


"""
Loss function for multi-exponential fit.
Returns the sum of squared differences between the data and the model.
"""
function mexp_loss(u, p)
    x = p[1]
    y = p[2]
    return norm((mexp(u, x) .- y) , p[3])
end

export expfit
"""
    expfit(
    input::ExperimentData{1},
    n::Union{Int, Vector{<:Real}}=1;
    solver= IPNewton(),
    normalize::Bool=true,
    L::Int = 2
)

Fit an n-exponential function to the data in the input struct.

The outut is an `ExpfitData` structure.

Arguments:

- `input` : 1D `ExperimentData`.
- `n` : Integer specifying the number of exponential terms.

Optional arguments:
- `solver` : Optim solver, default choice is IPNewton().
- `L` : An integer specifying which norm of the residuals you want to minimize (default is 2).

The `n` argument can also be a vector of initial parameter guesses, 
in which case it also determines the number of exponential terms used.
It should be of the form [a1, b1, a2, b2, ...], 
where a's are the amplitudes and b's are the parameters inside the exponentials.

The length of the vector in this case must be an even number.

The following examples might help to clarify: 

`expfit(input, [a,b])` -> mono-exponential fit with initial guess: a * exp.( (1/b) * x) 


`expfit(input, [a,b,c,d])` -> bi-exponential fit with initial guess: a * exp.( (1/b) * x) + c * exp.((1/d) * x) 


`expfit(input, [a,b,c,d,e,f])` -> tri-exponential fit with initial guess: a * exp.( (1/b) * x) + c * exp.((1/d) * x) + e * exp.((1/f) * x) 

(where a,b,c,d,e,f are numbers of your choice)

Numbers of parameters beyond tri-exponential can also be used, but it is not recommended.
"""
function expfit(
    input::ExperimentData{1},
    n::Union{Int, Vector{<:Real}}=1;
    solver= IPNewton(),
    scale::Bool = false,
    L::Int = 2,
)

    if scale 
        # copy so that original data is not mutated
        input = deepcopy(input)
        scale_to_one!(input.data)
    end

    x = input.axes[1]
    y = real.(input.data)

    if isa(n,Int) 

        # Make an educated guess of the initial parameters
        if x isa PFG
            u0 = [
                v for l in 1:n for v in ( 
                    abs(y[1])/(3*l -2) ,maximum(x)/100 * 10.0^(-(l/2 -0.5)) 
                ) 
            ]
        else 
            u0 = [
                v for l in 1:n for v in (
                    abs(y[1])/(3*l -2) ,maximum(x)/5 * 10.0^(-(l/2))
                )
            ]
        end

    elseif isa(n,Vector)

        u0 = Float64.(n)

    end

    u = optimize(
        u -> mexp_loss(u, (x, y, L)), 
        zeros(length(u0)), Inf .* ones(length(u0)), u0, 
        solver
    ).minimizer

    # Determine what's the x-axis of the seq (time or bfactor)
    x isa PFG ? x_ax = "b" : x_ax = "t"

    # Determine whether there should be a multiplication or a division in the exp() formula 
    x isa PFG ?  op = "*" : op = "/"

    un = sum([u[i] for i in 1:2:length(u)]) 

    # Get the fit's equation as a string
    eq = join(
        [
            (i == 1 ? "" : " + ") * 
            string(round(u[i], sigdigits=2)) * 
            " * exp(-$x_ax" * op * string(round(u[i+1], sigdigits=2)) * 
            ")" 
            for i in 1:2:length(u)
        ]
    )

    # Normalised version of the fit equation
    eqn = string( round(un, sigdigits=2)) * " * (" * 
        join(
            [
                (i == 1 ? "" : " + ") * 
                string(round(u[i] / un, sigdigits=2)) * 
                " * exp(-$x_ax" * op * string(round(u[i+1], sigdigits=2)) * ")" 
                for i in 1:2:length(u)
            ]
        ) * ")"

    if x isa IR
        unstr = string(round(un,sigdigits = 2))
        eq = unstr * " - 2 * (" * eq * ")" 
        eqn = unstr * " - 2 * " *  eqn   

    elseif x isa SR
        unstr = string(round(un,sigdigits = 2))
        eq = unstr * " -" * "(" * eq * ")" 
        eqn = unstr * " -" *  eqn   
    end

    return ExpfitData(
        input,
        u,
        u0,
        mexp(u, x) - y,
        eq,
        length(u0) == 2 ? eq : eqn,
        ""
    )
end
