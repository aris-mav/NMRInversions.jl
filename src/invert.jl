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
- `alpha` : can either be a single scalar value (e.g. 0.65) or an `AlphaOptimizer`.
- `solver` : is the `RegularizationSolver` used to solve the optimization problem.
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
    alpha::Union{Real,AlphaOptimizer}=GCV(),
    solver::Union{RegularizationSolver,Type{<:RegularizationSolver}}=BRD(),
    silent::Bool=true,
    scale::Bool=true,
) where {D}

    dims_to_invert = haskernel.(input.axes)

    if all(!, dims_to_invert)
        error("None of the dimensions are invertable.")
    end

    if scale
        # copy so that original data is not mutated
        input = deepcopy(input)
        scale_to_one!(input.data)
    end

    n_points = div(128, 2^(ndims(input) - 1))

    axes = ntuple(Val(ndims(input))) do i
        if dims_to_invert[i]
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

    # if dim should be inverted, select all of it
    # if not, iterate through its elements
    slices = Tuple(
        inv ? (Colon(),) : Base.axes(input, dim)
        for (dim, inv) in enumerate(dims_to_invert)
    )

    data = Array{Float64}(undef, length.(axes))
    alphas = Array{Float64}(undef, length.(axes))
    residuals = Array{Float64}(undef, size(input.data))

    for idx in Iterators.product(slices...)

        if !silent && !all(dims_to_invert)
            println("Inverting \
            $(nameof.(typeof.(input.axes))[findall(dims_to_invert)]), \
            at data$(_show_idx(idx))")
        end

        scale_kernel = count(dims_to_invert) < length(dims_to_invert) ?
                       true : false

        data[idx...], residuals[idx...], alphas[idx...] = _solve(
            input[idx...], axes[findall(dims_to_invert)],
            alpha, solver, silent, scale_kernel)
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
        residuals,
        alphas,
        ones(eltype(data), size(data)),
        [],
        "",
    )
end

function _solve(
    input::ExperimentData, axes,
    alpha::Union{Real,AlphaOptimizer},
    solver::RegularizationSolver,
    silent::Bool,
    kscale::Bool,
)

    D = ndims(input)

    if D > 2
        error("3D inversions not implemented yet, please submit an issue.")
    end

    k = create_kernel(input, axes, scale=kscale)

    if isa(alpha, Real)
        f, _ = solve_regularization(k.K, k.g, alpha, solver)
        α = alpha
    else
        f, _, α = find_alpha(k, solver, alpha; silent=silent)
    end

    Kf = if D == 1
        k.K_uncompressed[1] * f
    elseif D == 2
        k.K_uncompressed[1] * reshape(f, length.(axes)) * k.K_uncompressed[2]'
    end

    g = real.(input.data)

    r = Kf - g

    if isa(α, Real)
        alphas = fill(α, size(f))
    end

    return f, r, alphas
end

function _show_idx(idx)
    "[" * join((x isa Colon ? ":" : string(x) for x in idx), ",") * "]"
end
