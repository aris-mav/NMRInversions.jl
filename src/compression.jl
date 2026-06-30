# abstract type DataAxis{T} <: AbstractVector{T} end;

# struct ExperimentData{D}
#     axes::NTuple{D, DataAxis}
#     data::AbstractArray{<:Number, D}
# end

"""
    window_average(G::Matrix{T}, target_rows::Int) where T

Compress data in one of its dimensions by dividing into non-overlaping bins and 
averaging the values within. Returns the compressed data, as well as the 
precision matrix (inverse of the covariance matrix).
"""
function window_average(
    input::ExperimentData,
    target_length::Int,
    dims::Int = 1,
)
    x = input.axes[dims]
    data = input.data

    original_length = length(x)

    if original_length <= target_length
        error("Target length must be smaller than the current length.")
    end

    # Construct evenly spaced bin edges in axis coordinate space
    # These define how we partition the axis range into compression windows
    edges = range(first(x), last(x), length = target_length + 1)

    # Convert those axis-space edges into index-space boundaries
    # Each edge is mapped to the first index where x >= edge
    idx_edges = [searchsortedfirst(x, e) for e in edges]

    # Ensure the final boundary includes the full array
    idx_edges[end] = original_length + 1

    # Fix potential issues where repeated indices create empty bins
    # (can happen if x is non-uniform or has duplicates)
    for i in 2:length(idx_edges)
        idx_edges[i] = max(idx_edges[i], idx_edges[i-1])
    end

    # Compute the shape of the output array:
    # same as input, but with compressed length along `dims`
    new_size = Base.setindex(size(data), target_length, dims)

    # Create new axis container of same type as original axis
    # but with compressed length
    new_axis = typeof(x)(similar(x.x, target_length))

    # Preallocate output data array with correct shape
    new_data = similar(input.data, new_size)

    # Loop over each compression bin
    for b in 1:target_length

        # Define index range for this bin in original data
        lo = idx_edges[b]
        hi = idx_edges[b + 1] - 1

        # Safety check: ensure bin is not empty
        if lo > hi
            error("Encountered an empty bin. Try using fewer bins.")
        end

        # Convert boundaries into a Julia range for slicing
        inds = lo:hi

        # Compute new axis value as the mean coordinate of the original axis
        new_axis[b] = mean(view(x.x, inds))

        # Build a multidimensional index that selects this bin along `dims`
        # and selects all elements along other dimensions
        select = ntuple(d ->
            d == dims ? inds : Colon(),
            ndims(data)
        )

        # Map the current bin's target position into a multi-dimensional 
        # index slice, keeping the target dimension intact for uniform array shapes.
        out_select = ntuple(d ->
            d == dims ? (b:b) : Colon(),
            ndims(data)
        )

        # Compress the data within this bin's window by calculating the mean, 
        # then write the result into the corresponding slot of the downsampled array.
        new_data[out_select...] .= mean(view(data, select...), dims=dims)
    end

    # replace only the compressed axis, keep other axes unchanged
    return ExperimentData(
        Base.setindex(input.axes, new_axis, dims),
        new_data
    )
end
