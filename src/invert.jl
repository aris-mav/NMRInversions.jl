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
    axes::NTuple{
        D,Union{AbstractVector{<:Real},AbstractRange,Nothing}
    }=ntuple(i -> nothing, Val(D)),
    alpha::Union{Real,alpha_optimizer}=gcv(),
    solver::Union{
        regularization_solver,Type{<:regularization_solver}
    }=brd(),
    silent::Bool=false,
    scale::Bool=true,
) where {D}

    invertable_dims = haskernel.(input.axes)

    if all(!, invertable_dims)
        error("None of the dimensions are invertable.")
    end

    if scale
        # copy so that original data is not mutated
        input = deepcopy(input)
        scale_to_one!(input.data)
    end

    n_points = div(128, 2^(length(input.axes) - 1))

    axes = ntuple(Val(length(input.axes))) do i
        if invertable_dims[i]
            if isnothing(axes[i])
                x = input.axes[i].x
                lower, upper = x isa PFG ?
                               (minimum(x) * 7, maximum(x) * 10) :
                               (minimum(x) / 7, maximum(x) * 7)

                values = NMRInversions.logrange(lower, upper, n_points)
            else
                values = axes[i]
            end
            typeof(input.axes[i])(values)
        else
            input.axes[i]
        end
    end

    data, r, α =
        if all(invertable_dims)

            f, r, α = _solve(create_kernel(input, axes), alpha, solver, silent)
            data = reshape(f, length.(axes))

            data, r, α
        else
            # which dimensions get a fixed index each iteration 
            # (in slowest-to-fastest order)
            dims_to_iterate = reverse(findall(!, invertable_dims))

            # the actual index ranges for those dimensions, 
            # e.g. (1:1, 1:3) for dims (3,2)
            ranges = Base.axes.(Ref(input), dims_to_iterate)

            data = Array{Float64}(undef, length.(axes)...)

            for values in Iterators.product(ranges...)
                # `values` is a tuple like (i3, i2), 
                # one fixed index per iterated dim,
                # in the same order as dims_to_iterate

                # map each iterated dimension to its current fixed value
                fixed = Dict(dims_to_iterate .=> values)

                # build the full index: `:` where invertable, fixed[d] otherwise
                idx = ntuple(length(invertable_dims)) do d
                    invertable_dims[d] ? Colon() : fixed[d]
                end

                silent || println("Inverting slice $(_show_idx(idx))")

                data[idx...], r, α = _solve(
                    create_kernel(
                        input[idx...], axes[findall(invertable_dims)]
                    ),
                    alpha, solver, silent
                )
            end
            data, r, α
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
        data,
        r,
        α,
        ones(eltype(data), size(data)),
        [],
        "",
    )
end

function _solve(ker_struct::svd_kernel_struct,
    alpha::Real,
    solver::regularization_solver,
    _
)
    f, r = solve_regularization(
        ker_struct.K, ker_struct.g, alpha, solver
    )
    return f, r, alpha
end

function _solve(ker_struct::svd_kernel_struct,
    alpha::alpha_optimizer,
    solver::regularization_solver,
    silent
)
    find_alpha(ker_struct, solver, alpha; silent=silent)
end

function _show_idx(idx)
    "[" * join((x isa Colon ? ":" : string(x) for x in idx), ",") * "]"
end
