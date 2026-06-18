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

    for i in eachindex(X_mut)
        if isnothing(X_mut[i])
            x = axes[i].x
            if x isa PFG
                X_mut[i] = collect(NMRInversions.logrange(
                    minimum(x)*7 , maximum(x)*10 , 128
                )) 
            else
                X_mut[i] = collect(NMRInversions.logrange(
                    minimum(x)/7 , maximum(x)*7 , 128
                ))
            end
        end
    end

    X = (X_mut...,)

    ker_struct = create_kernel(axes, X, g)

    α = 0.0 #placeholder, will be replaced below 
    if alpha isa Real
        α = alpha
        f, r = solve_regularization(ker_struct.K, ker_struct.g, α, solver)
    else
        f, r, α = find_alpha(ker_struct, solver, alpha ; silent = silent)
    end

    if axes isa PFG
        X .= X ./ 1e9
    end

    for i in eachindex(X)
        if axes[i].x isa PFG
            X[i] .= X[i] ./ 1e9
        end
    end

    return InversionData{D}(
        axes, data, X, f, r, SNR, α, ones(eltype(data), size(data)), [], ""
    )

end

function invert(data::ExperimentData; kwargs...)
    return invert(data.axes, data.data; kwargs...)
end
