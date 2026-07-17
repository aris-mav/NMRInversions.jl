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
- `Spectrum`
- `FID` (Free induction decay)
- `FC` (Field cycling)
"""
abstract type DataAxis{T} <: AbstractVector{T} end;
export DataAxis

# Make DataAxis indexable
Base.size(E::DataAxis) = size(E.x)
Base.length(E::DataAxis) = length(E.x)
Base.IndexStyle(::Type{<:DataAxis}) = IndexLinear()
Base.getindex(E::DataAxis, i::Int) = E.x[i]
Base.setindex!(E::DataAxis, val, i::Int) = (E.x[i] = val)

function Base.getindex(E::T, I::AbstractArray{<:Integer}) where {T<:DataAxis}
    T(E.x[I])
end
function Base.getindex(E::T, I::AbstractRange{<:Integer}) where {T<:DataAxis}
    T(E.x[I])
end

# Promote integer inputs to float for any DataAxis subtype
function (::Type{T})(x::AbstractVector{<:Integer}) where {T<:DataAxis}
    T(float.(x))
end

# Define all the DataAxis types

"""
    IR{T<:Real} <: DataAxis{T}

`DataAxis` type for inversion recovery measurements.
The `x` field contains the corresponding time data.
"""
struct IR{T<:Real} <: DataAxis{T}
    x::AbstractVector{T}
end;
export IR

"""
    SR{T<:Real} <: DataAxis{T}

`DataAxis` type for saturation recovery measurements.
The `x` field contains the corresponding time data.
"""
struct SR{T<:Real} <: DataAxis{T}
    x::AbstractVector{T}
end;
export SR

"""
    CPMG{T<:Real} <: DataAxis{T}

`DataAxis` type for CPMG measurements.
The `x` field contains the corresponding time data.
"""
struct CPMG{T<:Real} <: DataAxis{T}
    x::AbstractVector{T}
end;
export CPMG

"""
    PFG{T<:Real} <: DataAxis{T}

`DataAxis` type for pulsed field gradient measurements.
The `x` field contains the corresponding b-factor data.
"""
struct PFG{T<:Real} <: DataAxis{T}
    x::AbstractVector{T}
end;
export PFG

"""
    Spectrum{T<:Real} <: DataAxis{T}

`DataAxis` type for NMR spectra.
The `x` field contains the corresponding ppm data.
"""
struct Spectrum{T<:Real} <: DataAxis{T}
    x::AbstractVector{T}
end;
export Spectrum

"""
    FID{T<:Real} <: DataAxis{T}

`DataAxis` type for free induction decay measurements.
The `x` field contains the corresponding time data.
"""
struct FID{T<:Real} <: DataAxis{T}
    x::AbstractVector{T}
end;
export FID

"""
    FC{T<:Real} <: DataAxis{T}

`DataAxis` type for field cycling (or "dispersion") measurements.
The `x` field contains the corresponding B0 data.
"""
struct FC{T<:Real} <: DataAxis{T}
    x::AbstractVector{T}
end;
export FC


# Define experiment data containing both Axes and data 
export ExperimentData
"""
    ExperimentData(axes, data)
Struct containing NMR data.
- `axes` should be a tuple of `DataAxis`s with the `x` data in each of them.
- `data` should be an aray of n-dimensions with the recorded data.
- `SNR` is the signal-to-noise ratio (calculated automatically for complex data)
- `W` is a tuple containing the precision matrix (inverse of the covariance 
matrix) in each of the dimensions. Relevant for compressed data, identity by 
default.
"""
struct ExperimentData{D}
    axes::NTuple{D,DataAxis}
    data::AbstractArray{<:Number,D}
    SNR::Real
    W::NTuple{D,Union{AbstractMatrix,UniformScaling}}
end

# Flexible outer constructor to catch dimension size mismatches
"""
    ExperimentData(
    axes::Tuple{Vararg{DataAxis}}, 
    data::AbstractArray{<:Number,D};
    compress = false
) where {D}

Constructor for `ExperimentalData`, calculating SNR and compressing the data if 
needed.
"""
function ExperimentData(
    axes::Tuple{Vararg{DataAxis}},
    data::AbstractArray{<:Number,D};
    compress=false
) where {D}

    if length(axes) != D
        throw(DimensionMismatch(
            """
            The number of axes (you gave $(length(axes))) \
            must match the dimensions of the data ($D).
            """
        ))
    end

    for i in 1:D
        if length(axes[i].x) != size(data, i)
            throw(DimensionMismatch(
                """
                The axis $(typeof(axes[i])) has length $(length(axes[i].x)), \
                but the `data` array has length $(size(data, i)) \
                at the corresponding dimension $i.
                """
            ))
        end
    end

    original_data = ExperimentData{D}(
        axes, data, calc_snr(axes, data), ntuple(i -> LinearAlgebra.I, D)
    )

    if compress
        return compress(original_data)
    else
        return original_data
    end
end

# Outer constructor for 1-dimensional case
function ExperimentData(axes::DataAxis, data::AbstractVector)
    return ExperimentData((axes,), data)
end

# Implement indexing
Base.length(d::ExperimentData) = length(d.data)
Base.size(d::ExperimentData) = size(d.data)
Base.size(d::ExperimentData, dim::Int) = size(d.data, dim)
Base.lastindex(d::ExperimentData) = lastindex(d.data)
Base.lastindex(d::ExperimentData, dim::Int) = size(d.data, dim)
Base.axes(d::ExperimentData, dim::Int) = axes(d.data, dim)
Base.ndims(d::ExperimentData) = ndims(d.data)

# Slice W matrix
_slice_W(W::AbstractMatrix, idx) = W[idx, idx]
# Or just return W if it is LinearAlgebra.I or similar
_slice_W(W::UniformScaling, _) = W

function Base.getindex(E::ExperimentData{D}, I...) where {D}

    idx = Base.to_indices(E.data, I)

    if length(idx) != D
        throw(ArgumentError(
            "Provide one index per dimension (expected $D, got $(length(idx))).")
        )
    end

    keep = findall(i -> !(idx[i] isa Integer), 1:D)
    new_D = length(keep)

    if new_D == 0
        return E.data[idx...] # fully scalar index -> return the raw value
    end

    new_axes = ntuple(j -> E.axes[keep[j]][idx[keep[j]]], new_D)
    new_W = ntuple(j -> _slice_W(E.W[keep[j]], idx[keep[j]]), new_D)

    return ExperimentData{new_D}(
        new_axes,
        E.data[idx...],
        E.SNR,
        new_W
    )
end

export InversionData
"""
    InversionData()

Output of the `invert` function, containing all the relevant information:

- `input` : `ExperimentData` structure that was used as input to the inversion function.
- `axes` : "x-axes" of the output (i.e., relaxation times, diffusion coefficients).
- `data` : the results of the inversion.
- `residuals` : an array containing the residuals of the fit.
- `alpha` : the smooting term used to produce the `data` array.
- `filter` : an array working as a mask to the `data` array (mostly for visualization).
- `selections` : polygon selections for the `data` array (mostly for visualization).
- `title` : a title descibing this object (mostly for visualization).
"""
mutable struct InversionData{D}
    input::ExperimentData{D}
    axes::NTuple{D,DataAxis}
    data::AbstractArray{<:Number,D}
    residuals::AbstractArray{<:Number,D}
    alphas::AbstractArray{<:Number,D}
    filter::AbstractArray{<:Number,D}
    selections::Vector{Vector{Vector}}
    title::String
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
