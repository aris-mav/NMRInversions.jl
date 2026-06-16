## The following are custom types, defined mainly for multiple dispatch purposes

"Superset of supported solvers for regularization."
abstract type regularization_solver end ; export regularization_solver

"Superset of supported methods to determine the α parameter in tikhonov
regularization."
abstract type alpha_optimizer end; export alpha_optimizer

"""
    NMRAxis{T} <: AbstractVector{T}

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
abstract type NMRAxis{T} <: AbstractVector{T} end

Base.size(E::NMRAxis) = size(E.x)
Base.getindex(E::NMRAxis, i::Int) = E.x[i]
Base.setindex!(E::NMRAxis, val, i::Int) = (E.x[i] = val)

"""
    IR{T<:Real} <: NMRAxis{T}

`NMRAxis` type for inversion recovery measurements.
The `x` field contains the corresponding time data.
"""
struct IR{T<:Real} <: NMRAxis{T} x::Vector{T} end; export IR

"""
    SR{T<:Real} <: NMRAxis{T}

`NMRAxis` type for saturation recovery measurements.
The `x` field contains the corresponding time data.
"""
struct SR{T<:Real} <: NMRAxis{T} x::Vector{T} end; export SR

"""
    CPMG{T<:Real} <: NMRAxis{T}

`NMRAxis` type for CPMG measurements.
The `x` field contains the corresponding time data.
"""
struct CPMG{T<:Real} <: NMRAxis{T} x::Vector{T} end; export CPMG

"""
    PFG{T<:Real} <: NMRAxis{T}

`NMRAxis` type for pulsed field gradient measurements.
The `x` field contains the corresponding b-factor data.
"""
struct PFG{T<:Real} <: NMRAxis{T} x::Vector{T} end; export PFG

"""
    Spectrum{T<:Real} <: NMRAxis{T}

`NMRAxis` type for NMR spectra.
The `x` field contains the corresponding ppm data.
"""
struct Spectrum{T<:Real} <: NMRAxis{T} x::Vector{T} end; export Spectrum

"""
    FID{T<:Real} <: NMRAxis{T}

`NMRAxis` type for free induction decay measurements.
The `x` field contains the corresponding time data.
"""
struct FID{T<:Real} <: NMRAxis{T} x::Vector{T} end; export FID

"""
    FC{T<:Real} <: NMRAxis{T}

`NMRAxis` type for field cycling (or "dispersion") measurements.
The `x` field contains the corresponding B0 data.
"""
struct FC{T<:Real} <: NMRAxis{T} x::Vector{T} end; export FC

"""
    NMRData(axes, data)
Struct containing NMR data.
- `axes` should be a tuple of `NMRAxis`s with the `x` data in each of them.
- `data` should be an aray of n-dimensions with the recorded data.
"""
struct NMRData{
    D,
    X <: NTuple{D, NMRAxis},
    Y <: AbstractArray{<:Number, D}
}
    axes::X
    data::Y
    # Strict inner constructor to enforce matching types
    function NMRData{D, X, Y}(axes::X, data::Y) where {D, X, Y}
        return new{D, X, Y}(axes, data)
    end
end

# Flexible outer constructor to catch dimension size mismatches
function NMRData(
    axes::Tuple{Vararg{NMRAxis}}, 
    data::AbstractArray{<:Number,D}
) where {D}
    len_x = length(axes)
    
    if len_x != D
        throw(DimensionMismatch(
            """
            The number of axes (you gave $len_x) \
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

    return NMRData{D, typeof(axes), typeof(data)}(axes, data)
end

# Outer constructor for 1-dimensional case
function NMRData(axes::NMRAxis, data::AbstractVector)
    return NMRData((axes,), data)
end

export InversionData
"""
    InversionData()

Output of the invert function for 1D pulse sequences.
A structure containing the following fields:
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PFG)
- `X`, the x values of the inversion results.
- `f`, the inversion results.
- `r`, the residuals.
- `SNR`, the signal-to-noise ratio.
- `α`, the regularization parameter.
- `selections` look like: `selections[sel_a[[x1,y1,..],[x2,y2,..],..], sel_b..]`
- `filter`, a vector representing the mask used to filter and scale the data.
- `title`, a title describing the data.

"""
mutable struct InversionData{
    D,
    X <: NTuple{D, NMRAxis},
    Y <: AbstractArray{<:Number, D}
}
    axes::X
    data::Y
    f::Y
    r::Y
    SNR::Real
    alpha::Real
    filter::Y
    selections::Vector{Vector{Vector}} 
    title::String
end

function InversionData(
    axes::Tuple{Vararg{NMRAxis}}, 
    data::AbstractArray{<:Number,D},
    f::AbstractArray{<:Number,D},
    r::AbstractArray{<:Number,D};
    SNR::Real = NaN,
    alpha::Real = NaN,
    filter::AbstractArray{<:Number,D} = ones(eltype(data), size(data)),
    selections::Vector{Vector{Vector{T}}} = Vector{Vector{Any}}[],
    title::String = "",
) where {D, T}
    
    for element in [data, f, r, filter]
        for i in 1:D
            if length(axes[i].x) != size(element, i)
                throw(DimensionMismatch(
                    """
                    The axis $i ($(typeof(axes[i]))) has a length of \
                    $(length(axes[i].x)), but the `data` array has length of \
                    $(size(data, i)) at the corresponding dimension $i.
                    """
                ))
            end
        end
    end

    return InversionData{D, typeof(axes), typeof(data)}(
        axes, data, f, r, SNR, alpha, filter, selections, title
    )
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
    X <: NTuple{D, NMRAxis},
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
