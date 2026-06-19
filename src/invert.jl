export invert
function invert(
    axes::NTuple{D, DataAxis},
    data::AbstractArray{<:Number, D};
    X::NTuple{ D, Any } = ntuple(i -> nothing, Val(D)),
    alpha::Union{Real, alpha_optimizer} = gcv(),
    solver::Union{regularization_solver, Type{<:regularization_solver}} = brd(),
    silent::Bool = false,
    scale::Bool = true,
) where {D}

    g = data

    if scale
        scale_to_one!(g)
    end

    isreal(data) ? SNR = NaN : SNR = calc_snr(data)

    X_mut = Union{Nothing, Vector{<:Real}}[X...]

    n_points = div(128, 2^(length(axes)-1) )

    for i in eachindex(X_mut)
        if isnothing(X_mut[i])
            x = axes[i].x
            if x isa PFG
                X_mut[i] = collect(NMRInversions.logrange(
                    minimum(x)*7 , maximum(x)*10 , n_points
                )) 
            else
                X_mut[i] = collect(NMRInversions.logrange(
                    minimum(x)/7 , maximum(x)*7 , n_points
                ))
            end
        end
    end

    X = (X_mut...,)

    ker_struct = create_kernel(axes, X, g)

    f, r, α = if alpha isa Real
        f, r = solve_regularization(ker_struct.K, ker_struct.g, alpha, solver)
        f, r, alpha
    else
        find_alpha(ker_struct, solver, alpha; silent = silent)
    end

    for i in eachindex(X)
        if axes[i] isa PFG
            # convert to SI units
            X[i] .= X[i] ./ 1e9
       end
    end

    return InversionData{D}(
        axes, data, X, 
        reshape(f, length.(X)), 
        r, SNR, α, 
        ones(eltype(data), size(data)), [], ""
    )

end

function invert(data::ExperimentData; kwargs...)
    return invert(data.axes, data.data; kwargs...)
end
