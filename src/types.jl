"Superset of supported solvers for regularization."
abstract type regularization_solver end ; export regularization_solver

"Superset of supported methods to determine the α parameter in tikhonov
regularization."
abstract type alpha_optimizer end; export alpha_optimizer

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
abstract type DataAxis{T} <: AbstractVector{T} end; export DataAxis

# Make DataAxis indexable
#
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

# Define all the DataAxis types

"""
    IR{T<:Real} <: DataAxis{T}

`DataAxis` type for inversion recovery measurements.
The `x` field contains the corresponding time data.
"""
struct IR{T<:Real} <: DataAxis{T} x::Vector{T} end; export IR

"""
    SR{T<:Real} <: DataAxis{T}

`DataAxis` type for saturation recovery measurements.
The `x` field contains the corresponding time data.
"""
struct SR{T<:Real} <: DataAxis{T} x::Vector{T} end; export SR

"""
    CPMG{T<:Real} <: DataAxis{T}

`DataAxis` type for CPMG measurements.
The `x` field contains the corresponding time data.
"""
struct CPMG{T<:Real} <: DataAxis{T} x::Vector{T} end; export CPMG

"""
    PFG{T<:Real} <: DataAxis{T}

`DataAxis` type for pulsed field gradient measurements.
The `x` field contains the corresponding b-factor data.
"""
struct PFG{T<:Real} <: DataAxis{T} x::Vector{T} end; export PFG

"""
    Spectrum{T<:Real} <: DataAxis{T}

`DataAxis` type for NMR spectra.
The `x` field contains the corresponding ppm data.
"""
struct Spectrum{T<:Real} <: DataAxis{T} x::Vector{T} end; export Spectrum

"""
    FID{T<:Real} <: DataAxis{T}

`DataAxis` type for free induction decay measurements.
The `x` field contains the corresponding time data.
"""
struct FID{T<:Real} <: DataAxis{T} x::Vector{T} end; export FID

"""
    FC{T<:Real} <: DataAxis{T}

`DataAxis` type for field cycling (or "dispersion") measurements.
The `x` field contains the corresponding B0 data.
"""
struct FC{T<:Real} <: DataAxis{T} x::Vector{T} end; export FC


# Define experiment data containing both Axis (x) and data (y)

"""
    ExperimentData(axes, data)
Struct containing NMR data.
- `axes` should be a tuple of `DataAxis`s with the `x` data in each of them.
- `data` should be an aray of n-dimensions with the recorded data.
"""
struct ExperimentData{D}
    axes::NTuple{D, DataAxis}
    data::AbstractArray{<:Number, D}
end

# Flexible outer constructor to catch dimension size mismatches
function ExperimentData(
    axes::Tuple{Vararg{DataAxis}}, 
    data::AbstractArray{<:Number,D}
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

    return ExperimentData{D}(axes, data)
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

function Base.getindex(E::ExperimentData{D}, I...) where {D}

    idx = Base.to_indices(E.data, I)
    
    return ExperimentData(
        ntuple(i -> E.axes[i][idx[i]], D),
        E.data[idx...]
    )
end

export InversionData
"""
    InversionData()

Output of the invert function for 1D pulse sequences.
A structure containing the following fields:
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PFG)
- `f`, the inversion results.
- `r`, the residuals.
- `SNR`, the signal-to-noise ratio.
- `α`, the regularization parameter.
- `filter`, a vector representing the mask used to filter/scale the data.
- `selections` look like: `selections[sel_a[[x1,y1,..],[x2,y2,..],..], sel_b..]`
- `title`, a title describing the data.

"""
mutable struct InversionData{D}
    axes::NTuple{D, DataAxis}
    data::AbstractArray{<:Number, D}
    X::NTuple{ D, Vector{Real} }
    f::AbstractArray{<:Number, D}
    r::AbstractArray
    SNR::Real
    alpha::Real
    filter::AbstractArray{<:Number, D}
    selections::Vector{Vector{Vector}} 
    title::String
end


export ExpfitData
"""
Output of the expfit function.
Structure containing information about multiexponential fits.

The fields are as follows:
- `u` : The fitted parameters for the `mexp` function.
- `u0` : The initial parameters for the `mexp` function.
- `r` : The residuals.
- `eq` : The equation of the fitted function.
- `eqn` : Same as above, normalised.
- `title` : A title describing the data.
"""
mutable struct ExpfitData{
    D,
    X <: NTuple{D, DataAxis},
    Y <: AbstractArray{<:Number, D}
}
    axes::X
    data::Y
    u::Vector
    u0::Vector
    r::Vector
    eq::String
    eqn::String
    title::String
end


# define structs for extension solvers
function solve_regularization(K, g, α, solver::regularization_solver)
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
struct jump_nnls <: regularization_solver 
    order::Int
    solver::Symbol
end
