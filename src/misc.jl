export scale_to_one!
"Scale data so that the max value is one."
function scale_to_one!(a::AbstractArray)
    a .= a ./ a[argmax(real(a))]
end

"Defining this here so that it works on julia 1.10"
function logrange(a,b,c)
    exp10.(range(log10(a),log10(b),c))
end

"""
Pass some immutable structure and return it with some field(s) changed.
"""
function change_field(x::T; kwargs...) where {T}
    # Convert keyword overrides to a dictionary for fast lookup
    kw = Dict(kwargs)

    # Get the fields in the correct positional order
    vals = [ get(kw, f, getfield(x, f)) for f in fieldnames(T) ]

    # Call the canonical positional constructor
    return T(vals...)
end

"""
Create a finite difference differentiation matrix of order `order` for a matrix of size `m x m`.
"""
function Γ(m::Int, order::Int)
    # Eilers, P. H. C. (2003). Analytical Chemistry, 75(14), 3631–3636. (Supporting Information)
    # Taken from RegularizationTools.jl
    if order == 0
        return Array{Float64}(LinearAlgebra.I, (m, m))
    end
    return diff(Γ(m, order - 1), dims=1)
end

"""
    calc_snr(data)
Calculate the Signal-to-Noise Ratio (SNR) from complex data,
where the real part is (mostly) signal and the imaginary part is (mostly) noise.
The STD of the latter half of the imaginary signal is used for the calculation 
(former half might contain signal residues as well). 
"""
function calc_snr(data::AbstractArray{<:Complex})

    real_data = real.(data)
    imag_data = imag.(data)

    half_index = floor(Int, size(imag_data, 1) / 2)
    end_index = size(imag_data, 1)

    noise = collect(selectdim(imag_data, 1, half_index:end_index))

    σ_n = sqrt(
        sum((noise .- sum(noise) / length(noise)) .^ 2) / (length(noise) - 1)
    )

    SNR = maximum(abs.(real_data)) / σ_n

    return SNR
end

# """
# Return a vector of matrices, containing the F for each selection polygon.
# """
# function selections(res::inv_out_2D)
#
#     dir = res.X_direct
#     indir = res.X_indirect
#     F = res.F
#
#     x = collect(1:length(indir))
#     y = collect(1:length(dir))
#
#     z = res.F' .* res.filter'
#
#     points = [[i, j] for i in x, j in y] #important, used by inpolygon later
#     mask = zeros(size(points))
#
#     polygons = res.selections 
#     selections = [zeros(size(F)) for _ in 1:length(polygons)]
#
#     for (i, polygon) in enumerate(polygons)
#         mask .= [PolygonOps.inpolygon(p, polygon; in=1, on=1, out=0) for p in points]
#         selections[i] .= mask .* z 
#     end
#
#     return selections
#
# end

"""
    autophase(data::ExperimentData; rotation::Real=0)

Correct the phase on experimental data.

The `rotation` optional argument should be an angle (rad),
and will be added to the automatic angle correction.
"""
function autophase(data::ExperimentData; rotation::Real=0)

    y = data.data

    # Declare which points should be the (near) maxima for each type 
    rules = [
        IR   => last,
        SR   => last,
        CPMG => first,
        PFG  => first,
    ]

    for (type, func) in rules

        # find the 1st dimension of the chosen type
        dim = findfirst(isa.(data.axes, type))
        dim == nothing && continue

        # index all elements (:) of the chosen dimension, 1 for other dimensions
        idx = ntuple(
            i -> i == dim ? (:) : 1, 
            ndims(y)
        )

        # choose the slice corresponding to the current dimension
        y_slice = y[idx...]

        # isolate n data points (not too many)
        n = min(length(y_slice) ÷ 5 , 10)
        y_slice = func(y_slice, n)

        # look at the angles of these points
        angles = angle.(y_slice)

        # get correction angle as the average of the above
        θ = sum(angles) / n

        # rotate data
        y .*= exp(-im * θ + rotation)

        # check if first or last point is the correct way around
        # (redundant but doesn't hurt)
        if real(func(y[idx...])) < 0
            y .= -y
        end

        return ExperimentData(data.axes, y)
    end

    throw("Autophase not implemented for $(typeof.(data.axes)).")

end

"""
    _selection_indices(r::InversionData{1})

Returns a vector of ranges which correspond to selections.

Can be used as: `r.data[NMRInversions._selection_indices(r)[1]]`
"""
function _selection_indices(r::InversionData{1})
    [
        begin
            _, low = findmin(x -> abs(x - s[1][1]), r.axes[1])
            _, high = findmin(x -> abs(x - s[2][1]), r.axes[1])
            low:high
        end 
        for s in r.selections
    ]
end

export weighted_averages
"""
    weighted_averages(r::InversionData{1} ; silent::Bool = false)
Return a vector with the weighted averages for 
the selections in the input structure, and a 
vector with the respective area fractions of 
these selections.
"""
function weighted_averages(r::InversionData{1} ; silent::Bool = false)

    distribution = r.data .* r.filter
    wa = Vector(undef, length(r.selections))
    areas = Vector(undef, length(r.selections))
    total_area = sum(distribution)

    for (i,idx) in enumerate(_selection_indices(r))

        # prevents total area above 100 if selections are touching
        idx = first(idx)+1:last(idx)

        area = sum(distribution[idx])
        wa[i] = distribution[idx]' * r.axes[1][idx] / area
        areas[i] = area / total_area

        lbl = if r.axes[1] isa Union{IR,SR}
            ["<T₁>","s"]
        elseif r.axes[1] isa CPMG
            ["<T₂>","s"]
        elseif r.axes[1] isa PFG
            ["<D>", "m²/s"]
        end

        if !silent
            display("Selection $(collect('a':'z')[i]) :")
            display(lbl[1] * " = $(round(wa[i], sigdigits=4)) "*lbl[2])
            display("Area = $(round(areas[i], sigdigits=4) * 100) %")
            display("")
        end
    end

    return wa, areas
end


"""
    weighted_averages(r::inv_out_2D)
Return two vectors with the weighted averages 
for the selections in the input structure, one for each dimension,
as well as a vector with the volume fractions of these selections.
"""
function weighted_averages(r::InversionData{2}; silent::Bool = false)

    wa_indir = Vector(undef, length(r.selections))
    wa_dir = Vector(undef, length(r.selections))
    volumes = Vector(undef, length(r.selections))

    z = r.data' .* r.filter'
    x = r.axes[2]
    y = r.axes[1]

    points = [[i, j] for i in x, j in y]
    mask = zeros(size(points))

    for (i,s) in enumerate(r.selections)

        mask .= [PolygonOps.inpolygon(p, s; in=1, on=1, out=0) for p in points]
        spo = mask .* z

        indir_dist = vec(sum(spo, dims=2))
        dir_dist = vec(sum(spo, dims=1))

        wa_indir[i] = indir_dist' *  y / sum(spo)
        wa_dir[i] = dir_dist' *  y / sum(spo)
        volumes[i] = sum(spo)/sum(z)

        labels = Dict(
            :IR => ["<T₁>","s"],
            :SR => ["<T₁>","s"],
            :CPMG => ["<T₂>","s"],
            :PFG => ["<D>", "m²/s"],
        )
        dir_lbl = labels[nameof(typeof(r.axes[1]))]
        ind_lbl = labels[nameof(typeof(r.axes[2]))]

        if !silent
            display("Selection $(collect('a':'z')[i]) :")
            display(ind_lbl[1] * " = $(round(wa_indir[i], sigdigits=4)) "*ind_lbl[2])
            display(dir_lbl[1] * " = $(round(wa_dir[i], sigdigits=4)) "*dir_lbl[2])
            if r.axes[2] isa Union{IR,SR} && r.axes[1] isa CPMG
                display("T₁/T₂ = $(round(wa_indir[i]/wa_dir[i], sigdigits=2)) ")
            end
            display("Volume = $(round(volumes[i], sigdigits=4) * 100) %")
            println()
        end
    end

    return wa_indir, wa_dir, volumes
end


export scale_filter!
"""
    scale_filter!(res::InversionData{1}, range)
Add selected range to the filter, scaling it to keep the integral of `f` constant.
"""
function scale_filter!(res::InversionData)
    integral = sum(res.data)
    new_integral = sum(res.data .* res.filter)
    scale = new_integral != 0 ? integral / new_integral : 0
    res.filter .= res.filter .* scale
end
