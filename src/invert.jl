export invert
"""
    invert(input::ExperimentData{D}; kwargs...)

Receives `ExperimentData`, returns `InversionData`.

The keyword arguments are:

- `axes` : a D-dimensional tuple where each element is an array containing the 
range of the output (e.g. relaxation times or diffusion coefficients) in 
every corresponding dimension. `nothing` can be passed instead of an array in 
any of the tuple elements, in which case a sensible array will be generated 
automatically (that is also the default option).
- `alpha` : can either be a single scalar value (e.g. 0.65) or an `alpha_optimizer`.
- `solver` : is the `regularization_solver` used to solve the optimization problem.
- `silent` : whether the function should be printing information while it runs 
(defaults to `false`).
- `scale` : whether the input data should be scaled to a maximum value of 1 
(you may also call that "normalizing" the data).
"""
function invert(
    input::ExperimentData{D};
    axes::NTuple{D, Union{AbstractVector{<:Real}, AbstractRange, Nothing}} = ntuple(i -> nothing, Val(D)),
    alpha::Union{Real, alpha_optimizer} = gcv(),
    solver::Union{regularization_solver, Type{<:regularization_solver}} = brd(),
    silent::Bool = false,
    scale::Bool = true,
) where {D}

    if scale 
        # copy so that original data is not mutated
        input = deepcopy(input)
        scale_to_one!(input.data)
    end

    g = input.data

    n_points = div(128, 2^(length(input.axes)-1) )

    axes = ntuple(Val(length(input.axes))) do i
        current_axis = axes[i]
        axis_type = typeof(input.axes[i])

        if isnothing(current_axis)
            x = input.axes[i].x
            lower, upper = x isa PFG ? 
            (minimum(x)*7, maximum(x)*10) : 
            (minimum(x)/7, maximum(x)*7)

            axis_type(collect(NMRInversions.logrange(lower, upper, n_points)))

        elseif current_axis isa AbstractRange
            axis_type(collect(current_axis))
        else
            axis_type(current_axis)
        end
    end

    ker_struct = create_kernel(input.axes, axes, g)

    f, r, α = if alpha isa Real
        f, r = solve_regularization(ker_struct.K, ker_struct.g, alpha, solver)
        f, r, alpha
    else
        find_alpha(ker_struct, solver, alpha; silent = silent)
    end

    for i in eachindex(axes)
        if input.axes[i] isa PFG
            # convert to SI units
            axes[i] .= axes[i] ./ 1e9
       end
    end

    data = reshape(f, length.(axes))

    return InversionData{D}(
        input, 
        axes, 
        data,
        r, 
        α, 
        calc_snr(input.data),
        ones(eltype(data), size(data)), 
        [], 
        "",
    )
end
