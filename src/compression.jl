export compress
"""
Compress data in one of its dimensions by dividing into non-overlaping bins and 
averaging the values within.
Returns the compressed data, as well as the precision matrix (inverse of the
covariance matrix) within a single `ExperimentData` structure.

- `input` : the `ExperimentData` you want to compress.

Keyword arguments:

- `dims` : which dimension you want to compress (defaults to the largest).
- `bins` : the final length of the `dims` (defaults to a power of 2).
- `scale` : output x-scale, choose between `:auto`, `:log`, `:linear`.
"""
function compress(
    input::ExperimentData;
    dims::Int=0,
    bins::Int=0,
    scale::Symbol=:auto,
)

    if dims == 0
        dims = argmax(length.(input.axes))
    end

    x = input.axes[dims]

    original_length = length(x)

    if bins <= 1
        # automatically select a sensible target value, max 64 points by default
        candidates = 2 .^ (1:6)
        bins = candidates[argmin(map(
            x -> x <= 0 ? Inf : x, # candidate can't be larger than original
            original_length .- candidates # find the closest candidate
        ))]
    end

    if original_length <= bins
        @info("Target length must be smaller than the current length. \
        Returning original data without compression.")
        return input
    end

    logscale = if scale == :auto
        x isa Union{IR,SR,CPMG,PFG} ? true : false
    elseif scale == :log
        true
    elseif scale == :linear
        false
    end

    return window_average(input, bins, dims, logscale)
end

function window_average(
    input::ExperimentData,
    nbins::Int,
    dims::Int,
    log::Bool,
)

    x = input.axes[dims]
    data = input.data
    original_length = length(x)

    # separate 1st point and the rest, so that the 1st point is preserve
    rest_x = view(x, 2:original_length)
    
    if log
        edges = NMRInversions.logrange(first(rest_x), last(rest_x), nbins)
    else
        edges = range(first(rest_x), last(rest_x), length=nbins)
    end

    # Map the remaining edges to index boundaries
    rest_idx_edges = [searchsortedfirst(x, e) for e in edges]
    
    # Combine the first point's boundary with the rest
    # Bin 1 starts at 1. Bin 2 starts at 2. 
    # Subsequent bins follow our computed edges.
    idx_edges = vcat([1, 2], rest_idx_edges[2:end])
    
    # Ensure the final boundary includes the full array
    idx_edges[end] = original_length + 1

    # Monotonicity check
    for i in 2:length(idx_edges)
        idx_edges[i] = max(idx_edges[i], idx_edges[i-1] + 1)
    end

    if idx_edges[end] > original_length + 1
        @warn """
        Original data size $original_length is too small for a target_length of
            $nbins. Capping boundaries."
        """
        for i in 1:length(idx_edges)
            if idx_edges[i] > original_length + 1
                idx_edges[i] = original_length + 1
            end
        end
    end

    bin_counts = [idx_edges[b+1] - idx_edges[b] for b in 1:nbins]

    Wdims = Diagonal(bin_counts)

    new_size = Base.setindex(size(data), nbins, dims)

    new_axis = typeof(x)(similar(x.x, nbins))
    new_data = similar(input.data, new_size)

    for b in 1:nbins

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
    new_W = Base.setindex(input.W, Wdims, dims)

    return ExperimentData(
        new_axes,
        new_data,
        input.SNR,
        new_W
    )
end
