export invert
"""
    invert(input::ExperimentData{D}; kwargs...)

Receives `ExperimentData`, returns `InversionData`.

The keyword arguments are:

- `dimensions` : a D-dimensional tuple where each element is an array containing the \
range of the output (e.g. relaxation times or diffusion coefficients) in \
every corresponding dimension. `nothing` can be passed instead of an array in \
any of the tuple elements, in which case a sensible array will be generated \
automatically (that is also the default option).
- `alpha` : can either be a single scalar value (e.g. 0.65) or an `AlphaOptimizer`.
- `solver` : is the `RegularizationSolver` used to solve the optimization problem.
- `silent` : whether the function should be printing information while it runs \
(defaults to `false`).
- `scale` : whether the input data should be scaled to a maximum value of 1 \
(you may also call that "normalizing" the data).
"""
function invert(
    input::ExperimentData{D};
    dimensions::NTuple{
        D,Union{AbstractVector{<:Real},AbstractRange,Nothing}
    }=ntuple(i -> nothing, Val(D)),
    alpha::Union{Real,AlphaOptimizer}=GCV(),
    solver::Union{RegularizationSolver,Type{<:RegularizationSolver}}=BRD(),
    silent::Bool=true,
    scale::Bool=true,
) where {D}

    dims_to_invert = haskernel.(input.dimensions)

    if all(!, dims_to_invert)
        error("None of the dimensions are invertable.")
    end

    if scale
        # copy so that original data is not mutated
        input = deepcopy(input)
        scale_to_one!(input.signal)
    end

    n_points = div(128, 2^(ndims(input) - 1))

    dimensions = ntuple(Val(ndims(input))) do i
        if dims_to_invert[i]
            if isnothing(dimensions[i])
                x = input.dimensions[i].x
                lower, upper = x isa PFG ?
                               (minimum(x) * 7, maximum(x) * 10) :
                               (minimum(x) / 7, maximum(x) * 7)

                values = NMRInversions.logrange(lower, upper, n_points)
            else
                values = dimensions[i]
            end
            Base.typename(typeof(input.dimensions[i])).wrapper(values)
        else
            input.dimensions[i]
        end
    end

    # if dim should be inverted, select all of it
    # if not, iterate through its elements
    slices = Tuple(
        inv ? (Colon(),) : Base.axes(input, dim)
        for (dim, inv) in enumerate(dims_to_invert)
    )

    data = Array{Float64}(undef, length.(dimensions))
    alphas = Array{Float64}(undef, length.(dimensions))
    residuals = Array{Float64}(undef, size(input.signal))

    for idx in Iterators.product(slices...)

        if !silent && !all(dims_to_invert)
            println("Inverting \
            $(nameof.(typeof.(input.dimensions))[findall(dims_to_invert)]), \
            at data$(_show_idx(idx))")
        end

        data[idx...], residuals[idx...], alphas[idx...] = _solve(
            input[idx...], dimensions[findall(dims_to_invert)], alpha, solver, silent)
    end

    for i in eachindex(dimensions)
        if input.dimensions[i] isa PFG
            # convert to SI units
            dimensions[i] .= dimensions[i] ./ 1e9
        end
    end

    return InversionData{D}(
        input,
        dimensions,
        data,
        residuals,
        alphas,
        ones(eltype(data), size(data)),
        [],
        "",
    )
end

function _solve(
    input::ExperimentData, dimensions,
    alpha::Union{Real,AlphaOptimizer},
    solver::RegularizationSolver,
    silent::Bool,
)

    D = ndims(input)

    if D > 2
        error("3D inversions not implemented yet, please submit an issue.")
    end

    k = create_kernel(input, dimensions)

    if isa(alpha, Real)
        f, _ = solve_regularization(k.K, k.g, alpha, solver)
        α = alpha
    else
        f, _, α = find_alpha(k, solver, alpha; silent=silent)
    end

    Kf = if D == 1
        k.K_uncompressed[1] * f
    elseif D == 2
        k.K_uncompressed[1] * reshape(f, length.(dimensions)) * k.K_uncompressed[2]'
    end

    g = real.(input.signal)

    r = Kf - g

    if isa(α, Real)
        alphas = fill(α, size(f))
    end

    return f, r, alphas
end

function _show_idx(idx)
    "[" * join((x isa Colon ? ":" : string(x) for x in idx), ",") * "]"
end
