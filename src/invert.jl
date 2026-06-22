export invert
function invert(
    input::ExperimentData{D};
    axes::NTuple{D, Any} = ntuple(i -> nothing, Val(D)),
    alpha::Union{Real, alpha_optimizer} = gcv(),
    solver::Union{regularization_solver, Type{<:regularization_solver}} = brd(),
    silent::Bool = false,
    scale::Bool = true,
) where {D}

    g = input.data

    if scale
        scale_to_one!(g)
    end


    axes_mut = Union{Nothing, AbstractVector{<:Real}}[axes...]

    n_points = div(128, 2^(length(input.axes)-1) )

    for i in eachindex(axes_mut)
        if isnothing(axes_mut[i])
            x = input.axes[i].x
            lower, upper = if x isa PFG
                minimum(x)*7 , maximum(x)*10
            else
                minimum(x)/7 , maximum(x)*7
            end
        end
        axes_mut[i] = typeof(input.axes[i])(
            collect( NMRInversions.logrange(lower, upper, n_points) )
        )
    end

    axes = (axes_mut...,)

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

    return InversionData{D}(
        input, 
        axes, 
        reshape(f, length.(axes)), 
        r, 
        α, 
        isreal(input.data) ? NaN : calc_snr(input.data),
        ones(eltype(f), size(f)), 
        [], 
        "",
    )

end
