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

    # Bin edges in axis coordinates
    edges = range(first(x), last(x), length = target_length + 1)

    idx_edges = [searchsortedfirst(x, e) for e in edges]
    idx_edges[end] = original_length + 1

    # Prevent empty bins if duplicate indices occur
    for i in 2:length(idx_edges)
        idx_edges[i] = max(idx_edges[i], idx_edges[i-1])
    end

    new_axis_values = similar(x.x, target_length)

    new_size = collect(size(data))
    new_size[dims] = target_length

    new_data = Array{eltype(data)}(undef, Tuple(new_size))

    for b in 1:target_length
        lo = idx_edges[b]
        hi = idx_edges[b + 1] - 1

        if lo > hi
            error("Encountered an empty bin. Try using fewer bins.")
        end

        inds = lo:hi

        new_axis_values[b] = mean(view(x.x, inds))

        select = ntuple(d ->
            d == dims ? inds : Colon(),
            ndims(data)
        )

        out_select = ntuple(d ->
            d == dims ? b : Colon(),
            ndims(data)
        )

        new_data[out_select...] .= dropdims(
        mean(view(data, select...), dims=dims),
        dims=dims,
        )
    end
    new_axis = typeof(x)(new_axis_values)
    new_axes = Base.setindex(input.axes, new_axis, dims)

    return ExperimentData(new_axes, new_data)
end
