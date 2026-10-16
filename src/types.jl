"Superset of supported solvers for regularization."
abstract type RegularizationSolver end;
export RegularizationSolver

"Superset of supported methods to determine the α parameter in tikhonov
regularization."
abstract type AlphaOptimizer end;
export AlphaOptimizer

"""
    DataAxis{T} <: AbstractVector{T}

Superset for abstract NMR experiment types.
Subtypes are:
- `IR` (Inversion recovery)
- `SR` (Saturation recovery)
- `CPMG`
- `PFG` (Pulsed field gradient)
- `T1`
- `T2`
- `D`
- `Spectrum`
- `FID` (Free induction decay)
- `FC` (Field cycling)
"""
abstract type DataAxis{T} <: AbstractVector{T} end
export DataAxis

# define all the DataAxis types
for (S, desc) in (
    (:IR, "inversion recovery measurements. `x` contains time data (in seconds)."),
    (:SR, "saturation recovery measurements. `x` contains time data (in seconds)."),
    (:CPMG, "CPMG measurements. `x` contains time data (in seconds)."),
    (:PFG, "pulsed field gradient measurements. `x` contains b-factor data."),
    (:T1, "T1 relaxation time distributions. `x` contains time data (in seconds)."),
    (:T2, "T2 relaxation time distributions. `x` contains time data (in seconds)."),
    (:D, "distributions of diffusion coefficients. `x` contains diffusion coefficients."),
    (:Spectrum, "NMR spectra. `x` contains ppm data."),
    (:FID, "free induction decay measurements. `x` contains time data (in seconds)."),
    (:FC, "field cycling (\"dispersion\") measurements. `x` contains B0 data."),
)
    @eval begin
        @doc "    $($(string(S))){T<:Real} <: DataAxis{T}\n\n`DataAxis` type for " * $desc
        struct $S{T<:Real} <: DataAxis{T}
            x::Vector{T}
            # accept any AbstractVector (ranges, views, ...);
            $S{T}(x::AbstractVector) where {T<:Real} = new{T}(convert(Vector{T}, x))
        end
        # promote integers to floats
        $S(x::AbstractVector{T}) where {T<:Real} = $S{float(T)}(x)
        export $S
    end
end

# make DataAxis indexable
Base.size(E::DataAxis) = size(E.x)
Base.length(E::DataAxis) = length(E.x)
Base.IndexStyle(::Type{<:DataAxis}) = IndexLinear()
Base.getindex(E::DataAxis, i::Int) = E.x[i]
Base.setindex!(E::DataAxis, val, i::Int) = (E.x[i] = val)

Base.getindex(E::DataAxis, ::Base.Slice) = E   # `:` keeps the axis as is

function Base.getindex(E::DataAxis, I::AbstractArray{<:Integer})
    typeof(E).name.wrapper(E.x[I])
end


"""
    invert_axis(sampling_axis, values)
"""
invert_axis(::CPMG, v) = T2(v)
invert_axis(::Union{IR,SR}, v) = T1(v)
invert_axis(::PFG, v) = D(v)
invert_axis(::A, v) where {A<:DataAxis} = throw(ArgumentError(
    "no inversion defined for $(nameof(A)) axes"))
# Define experiment data containing both Axes and data 
export ExperimentData
"""
    ExperimentData(dimensions, signal)

Struct containing NMR data in its fields:

- `dimensions` should be a tuple of `DataAxis`s with the `x` data in each of them.
- `signal` should be an aray of n-dimensions with the recorded signal.
- `SNR` is the signal-to-noise ratio (calculated automatically for complex signals)
- `W` is a tuple containing the precision matrix (inverse of the covariance 
matrix) in each of the dimensions. Relevant for compressed data, identity by 
default.
"""
mutable struct ExperimentData{D}
    dimensions::NTuple{D,DataAxis}
    signal::AbstractArray{<:Number,D}
    SNR::Real
    W::NTuple{D,Union{AbstractMatrix,UniformScaling}}
end

# Flexible outer constructor to catch dimension size mismatches
"""
    ExperimentData(
    dimensions::Tuple{Vararg{DataAxis}}, 
    signal::AbstractArray{<:Number,D};
    compress = false
) where {D}

Constructor for `ExperimentalData`, calculating SNR and compressing the data if 
needed.
"""
function ExperimentData(
    dimensions::Tuple{Vararg{DataAxis}},
    signal::AbstractArray{<:Number,D};
    compress=false
) where {D}

    if length(dimensions) != D
        throw(DimensionMismatch(
            """
            The number of dimensions (you gave $(length(dimensions))) \
            must match the dimensions of the signal ($D).
            """
        ))
    end

    for i in 1:D
        if length(dimensions[i].x) != size(signal, i)
            throw(DimensionMismatch(
                """
                The axis $(typeof(dimensions[i])) has length $(length(dimensions[i].x)), \
                but the `signal` array has length $(size(signal, i)) \
                at the corresponding dimension $i.
                """
            ))
        end
    end

    original_data = ExperimentData{D}(
        dimensions, signal, calc_snr(dimensions, signal), ntuple(i -> LinearAlgebra.I, D)
    )

    if compress
        return compress(original_data)
    else
        return original_data
    end
end

# Outer constructor for 1-dimensional case
function ExperimentData(dimensions::DataAxis, signal::AbstractVector)
    return ExperimentData((dimensions,), signal)
end

# Implement indexing
Base.length(d::ExperimentData) = length(d.signal)
Base.size(d::ExperimentData) = size(d.signal)
Base.size(d::ExperimentData, dim::Int) = size(d.signal, dim)
Base.lastindex(d::ExperimentData) = lastindex(d.signal)
Base.lastindex(d::ExperimentData, dim::Int) = size(d.signal, dim)
Base.axes(d::ExperimentData, dim::Int) = axes(d.signal, dim)
Base.ndims(d::ExperimentData) = ndims(d.signal)

# Slice W matrix
_slice_W(W::AbstractMatrix, idx) = W[idx, idx]
# Or just return W if it is LinearAlgebra.I or similar
_slice_W(W::UniformScaling, _) = W

function Base.getindex(E::ExperimentData{D}, I...) where {D}

    idx = Base.to_indices(E.signal, I)

    if length(idx) != D
        throw(ArgumentError(
            "Provide one index per dimension (expected $D, got $(length(idx))).")
        )
    end

    keep = findall(i -> !(idx[i] isa Integer), 1:D)
    new_D = length(keep)

    if new_D == 0
        return E.signal[idx...] # fully scalar index -> return the raw value
    end

    new_dimensions = ntuple(j -> E.dimensions[keep[j]][idx[keep[j]]], new_D)
    new_W = ntuple(j -> _slice_W(E.W[keep[j]], idx[keep[j]]), new_D)

    return ExperimentData{new_D}(
        new_dimensions,
        E.signal[idx...],
        E.SNR,
        new_W
    )
end

function Base.show(io::IO, ::MIME"text/plain", d::ExperimentData{N}) where N
    print(io, "ExperimentData ($(join(size(d), "×")) points)")
end

export InversionData
"""
    InversionData()

Output of the `invert` function, containing all the relevant information in its
    fields:

- `input` : `ExperimentData` structure that was used as input to the inversion function.
- `dimensions` : "x-axes" of the output (i.e., relaxation times, diffusion coefficients).
- `distribution` : the results of the inversion.
- `residuals` : an array containing the residuals of the fit.
- `alpha` : the smooting term used to produce the `distribution` array.
- `filter` : an array working as a mask to the `distribution` array (mostly for visualization).
- `selections` : polygon selections for the `distribution` array (mostly for visualization).
- `title` : a title descibing this object (mostly for visualization).
"""
mutable struct InversionData{D}
    input::ExperimentData{D}
    dimensions::NTuple{D,DataAxis}
    distribution::AbstractArray{<:Number,D}
    residuals::AbstractArray{<:Number,D}
    alphas::AbstractArray{<:Number,D}
    filter::AbstractArray{<:Number,D}
    selections::Vector{Vector{Vector}}
    title::String
end

# Implement indexing
Base.length(d::InversionData) = length(d.distribution)
Base.size(d::InversionData) = size(d.distribution)
Base.size(d::InversionData, dim::Int) = size(d.distribution, dim)
Base.lastindex(d::InversionData) = lastindex(d.distribution)
Base.lastindex(d::InversionData, dim::Int) = size(d.distribution, dim)
Base.axes(d::InversionData, dim::Int) = axes(d.distribution, dim)
Base.ndims(d::InversionData) = ndims(d.distribution)

function Base.show(io::IO, ::MIME"text/plain", d::InversionData{N}) where N
    print(io, "InversionData ($(join(size(d), "×")) points)")
end

export ExpfitData
"""
Output of the `expfit` function.
Structure containing information about multiexponential fits.

The fields are as follows:
- `u` : The fitted parameters for the `mexp` function.
- `u0` : The initial parameters for the `mexp` function.
- `residuals` : The residuals of the fit.
- `eq` : The equation of the fitted function, printed as a string.
- `eqn` : Same as above, normalised.
- `title` : A title describing the data.
"""
mutable struct ExpfitData
    input::ExperimentData{1}
    u::AbstractVector
    u0::AbstractVector
    residuals::AbstractVector
    eq::String
    eqn::String
    title::String
end


# define structs for extension solvers
function solve_regularization(K, g, α, solver::RegularizationSolver)
    error("The package needed for $(typeof(solver)) is not loaded.")
end


export jump_nnls
"""
    jump_nnls(L, solver)
Jump non-negative least squares method for tikhonov (L2) regularization, 
implemented using the JuMP package extension.
All around effective, but can be slow for large problems, such as 2D inversions, 
unless a powerful solver like gurobi is used.

- `L` determines the tikhonov matrix.\
Can be either a matrix or an integer `n`. The latter will create a finite difference matrix\
of order `n`, which means that the penalty term will be the `n`'th derivative of the distribution.\
Defaults to `0`, where `L` becomes the Identity matrix.
- `solver` is passed directly into JuMP, see its documentation for available options. 

This one is still experimental and not well-tested, 
please submit an issue if you encounter any difficulties.
"""
struct jump_nnls <: RegularizationSolver
    order::Int
    solver::Symbol
end
