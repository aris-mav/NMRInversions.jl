export window_average
"""
    window_average(
    input::ExperimentData,
    target_length::Int=128,
    dims::Int=0;
    log::Bool=false,
)

Compress data in one of its dimensions by dividing into non-overlaping bins and 
averaging the values within.
Returns the compressed data, as well as the precision matrix (inverse of the
covariance matrix) within a single `ExperimentData` structure.
If dims is not provided, the largest dimension is chosen by default.
"""
function window_average(
    input::ExperimentData,
    target_length::Int=0,
    dims::Int=0;
    log::Bool=true,
)

    if dims == 0
        dims = argmax(length.(input.axes))
    end

    x = input.axes[dims]
    data = input.data

    original_length = length(x)

    if target_length <= 1
        # automatically select a sensible target value
        candidates = 2 .^ (1:10)
        target_length = candidates[argmin(map(
            x -> x < 0 ? Inf : x, # candidate can't be larger than original
            original_length .- candidates # find the closest candidate
        ))]
    end

    if original_length <= target_length
        @info("Target length must be smaller than the current length. \
        Returning original data without compression.")
        return input
    end

    # Construct evenly spaced bin edges in axis coordinate space
    # These define how we partition the axis range into compression windows
    if log
        edges = NMRInversions.logrange(first(x), last(x), target_length + 1)
    else
        edges = range(first(x), last(x), length=target_length + 1)
    end

    # Convert those axis-space edges into index-space boundaries
    # Each edge is mapped to the first index where x >= edge
    idx_edges = [searchsortedfirst(x, e) for e in edges]

    # Ensure the final boundary includes the full array
    idx_edges[end] = original_length + 1

    for i in 2:length(idx_edges)
        idx_edges[i] = max(idx_edges[i], idx_edges[i-1] + 1)
    end

    if idx_edges[end] > original_length + 1
        @warn """
        Original data size $original_length is too small for a target_length of
            $target_length. Capping boundaries."
        """
        for i in 1:length(idx_edges)
            if idx_edges[i] > original_length + 1
                idx_edges[i] = original_length + 1
            end
        end
    end

    bin_counts = [idx_edges[b+1] - idx_edges[b] for b in 1:target_length]
    W = Diagonal(bin_counts)

    new_size = Base.setindex(size(data), target_length, dims)
    new_axis = typeof(x)(similar(x.x, target_length))
    new_data = similar(input.data, new_size)

    for b in 1:target_length

        # Define index range for this bin in original data
        lo = idx_edges[b]
        hi = idx_edges[b+1] - 1

        # Safety check: ensure bin is not empty
        if lo > hi
            @warn("Encountered an empty bin. Returning original input.")
            return input
        end

        # Convert boundaries into a Julia range for slicing
        idx = lo:hi

        # Compute new axis value as the mean coordinate of the original axis
        new_axis[b] = mean(view(x.x, idx))

        # Build a multidimensional index that selects this bin along `dims`
        # and selects all elements along other dimensions
        select = ntuple(d ->
                d == dims ? idx : Colon(),
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

    new_axes = Base.setindex(input.axes, new_axis, dims)

    return ExperimentData(
        new_axes,
        new_data,
        calc_snr(new_axes, new_data),
        Base.setindex(input.W, W, dims),
    )

end
