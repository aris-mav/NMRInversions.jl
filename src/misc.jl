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
function calc_snr(data::AbstractMatrix{<:Complex}) # For matrix data

    real_data = real.(data)
    imag_data = imag.(data)

    noise = imag_data[(floor(Int64, (size(imag_data, 1) / 2))):end, :]
    σ_n = sqrt(sum((noise .- sum(noise) / length(noise)) .^ 2) / (length(noise) - 1))
    SNR = maximum(abs.(real_data)) / σ_n

    return SNR
end
function calc_snr(data::AbstractVector{<:Complex}) # For vector data

    real_data = real.(data)
    imag_data = imag.(data)

    noise = imag_data[(floor(Int64, (size(imag_data, 1) / 2))):end]
    σ_n = sqrt(sum((noise .- sum(noise) / length(noise)) .^ 2) / (length(noise) - 1))
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

function angle_correction(y::AbstractVecOrMat)
    n = min(floor(Int, length(y) / 5) , 15)
    θ = sum(angle.(y[1:n])) / n
    return θ
end

# function autophase(data::input1D; rotation::Real=0)
#
#     y_phased = data.y .* exp(-im * (angle_correction(data.y)+rotation) )
#
#     if rotation == 0
#
#         d = real.(y_phased[1] - y_phased[end])
#
#         if data.seq in (CPMG, PFG) && d < 0
#             return autophase(data, rotation = pi)
#         elseif data.seq in (IR, SR) && d > 0
#             return autophase(data, rotation = pi)
#         end
#     end
#
#     return input1D(data.seq, data.x, y_phased)
# end
#
# function autophase(data::input2D; rotation::Real=0)
#
#     y_phased = data.data .* exp(-im * (angle_correction(data.data)+rotation) )
#
#     if rotation == 0
#
#         d = real.(y_phased[1,1] - y_phased[1,end])
#
#         if data.seq in (PFGCPMG, CPMGCPMG) && d < 0
#             return autophase(data, rotation = pi)
#         elseif data.seq in (IRCPMG, SRCPMG) && d > 0
#             return autophase(data, rotation = pi)
#         end
#     end
#
#     return input2D(data.seq, data.x_direct, data.x_indirect, y_phased)
# end


# export weighted_averages
# """
#     weighted_averages(r::inv_out_1D)
# Return a vector with the weighted averages for 
# the selections in the input structure, and a 
# vector with the respective area fractions of 
# these selections.
# """
# function weighted_averages(r::inv_out_1D ; silent::Bool = false)
#
#     distribution = r.f .* r.filter
#     wa = Vector(undef, length(r.selections))
#     areas = Vector(undef, length(r.selections))
#     total_area = sum(distribution)
#
#     for (i,s) in enumerate(r.selections)
#         area = sum(distribution[s[1]:s[2]-1])
#         wa[i] = distribution[s[1]:s[2]-1]' * r.X[s[1]:s[2]-1] / area
#         areas[i] = area / total_area
#
#         lbl = if r.seq in [IR, SR]
#             ["<T₁>","s"]
#         elseif r.seq == CPMG
#             ["<T₂>","s"]
#         elseif r.seq == PFG
#             ["<D>", "m²/s"]
#         end
#
#         if !silent
#             display("Selection $(collect('a':'z')[i]) :")
#             display(lbl[1] * " = $(round(wa[i], sigdigits=4)) "*lbl[2])
#             display("Area = $(round(areas[i], sigdigits=4) * 100) %")
#             display("")
#         end
#     end
#
#     return wa, areas
# end
#
#
# """
#     weighted_averages(r::inv_out_2D)
# Return two vectors with the weighted averages 
# for the selections in the input structure, one for each dimension,
# as well as a vector with the volume fractions of these selections.
# """
# function weighted_averages(r::inv_out_2D; silent::Bool = false)
#
#     wa_indir = Vector(undef, length(r.selections))
#     wa_dir = Vector(undef, length(r.selections))
#     volumes = Vector(undef, length(r.selections))
#
#     z = r.F' .* r.filter'
#     x = r.X_indirect
#     y = r.X_direct
#
#     points = [[i, j] for i in x, j in y]
#     mask = zeros(size(points))
#
#     for (i,s) in enumerate(r.selections)
#
#         mask .= [PolygonOps.inpolygon(p, s; in=1, on=1, out=0) for p in points]
#         spo = mask .* z
#
#         indir_dist = vec(sum(spo, dims=2))
#         dir_dist = vec(sum(spo, dims=1))
#
#         wa_indir[i] = indir_dist' *  r.X_indirect / sum(spo)
#         wa_dir[i] = dir_dist' *  r.X_direct / sum(spo)
#         volumes[i] = sum(spo)/sum(z)
#
#         dir_lbl = ["<T₂>","s"]
#         ind_lbl = if r.seq in (IRCPMG, SRCPMG)
#             ["<T₁>","s"]
#         elseif r.seq == PFGCPMG
#             ["<D>", "m²/s"]
#         end
#
#         if !silent
#             display("Selection $(collect('a':'z')[i]) :")
#             display(ind_lbl[1] * " = $(round(wa_indir[i], sigdigits=4)) "*ind_lbl[2])
#             display(dir_lbl[1] * " = $(round(wa_dir[i], sigdigits=4)) "*dir_lbl[2])
#             if r.seq in (IRCPMG, SRCPMG)
#                 display("T₁/T₂ = $(round(wa_indir[i]/wa_dir[i], sigdigits=2)) ")
#             end
#             display("Volume = $(round(volumes[i], sigdigits=4) * 100) %")
#             println()
#         end
#     end
#
#     return wa_indir, wa_dir, volumes
# end


# export filter_range!
# """
#     filter_range!(res::inv_res_1D , range)
# Apply selected range to the filter, scaling it to keep the integral of `f` constant.
# """
# function filter_range!(res::inv_out_1D, range)
#
#     integral = sum(res.f)
#     res.filter[range[1]:range[2]] .= 0
#     new_integral = sum(res.f .* res.filter)
#     scale = new_integral != 0 ? integral / new_integral : 0
#     res.filter .= res.filter .* scale
#
# end



