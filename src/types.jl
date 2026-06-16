## The following are custom types, defined mainly for multiple dispatch purposes

"Superset of supported solvers for regularization."
abstract type regularization_solver end ; export regularization_solver

"Superset of supported methods to determine the α parameter in tikhonov
regularization."
abstract type alpha_optimizer end; export alpha_optimizer

"Superset of abstract NMR experiment types."
abstract type NMRExperiment end

# The reason why we don't do IR(x,y) is to be able to compose multidimensional
# experiments where x axes are independent but y data correlated
"""
    IR{T<:Real} <: NMRExperiment

`NMRExperiment` type for inversion recovery measurements.
The `x` field contains the corresponding time data.
"""
struct IR{T<:Real} <: NMRExperiment x::Vector{T} end; export IR

"""
    SR{T<:Real} <: NMRExperiment

`NMRExperiment` type for saturation recovery measurements.
The `x` field contains the corresponding time data.
"""
struct SR{T<:Real} <: NMRExperiment x::Vector{T} end; export SR

"""
    CPMG{T<:Real} <: NMRExperiment

`NMRExperiment` type for CPMG measurements.
The `x` field contains the corresponding time data.
"""
struct CPMG{T<:Real} <: NMRExperiment x::Vector{T} end; export CPMG

"""
    PFG{T<:Real} <: NMRExperiment

`NMRExperiment` type for pulsed field gradient measurements.
The `x` field contains the corresponding b-factor data.
"""
struct PFG{T<:Real} <: NMRExperiment x::Vector{T} end; export PFG

"""
    Spectrum{T<:Real} <: NMRExperiment

`NMRExperiment` type for NMR spectra.
The `x` field contains the corresponding ppm data.
"""
struct Spectrum{T<:Real} <: NMRExperiment x::Vector{T} end; export Spectrum

"""
    FID{T<:Real} <: NMRExperiment

`NMRExperiment` type for free induction decay measurements.
The `x` field contains the corresponding time data.
"""
struct FID{T<:Real} <: NMRExperiment x::Vector{T} end; export FID

"""
    Dispersion{T<:Real} <: NMRExperiment

`NMRExperiment` type for field cycling dispersion measurements.
The `x` field contains the corresponding B0 data.
"""
struct Dispersion{T<:Real} <: NMRExperiment x::Vector{T} end; export Dispersion


"""
Struct containing NMR data.
- `axes` should be a tuple of `NMRExperiment`s with the `x` data in each of them.
- `data` should be an aray of n-dimensions with the recorded data.
"""
struct NMRData{
    D,
    X <: NTuple{D, NMRExperiment},
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
    axes::Tuple{Vararg{NMRExperiment}}, 
    data::AbstractArray{<:Number,D}
) where {D}
    len_x = length(axes)
    
    if len_x != D
        throw(DimensionMismatch(
            """
            The number of experiments (you gave $len_x) must match 
            the number of array dimensions (you gave $D).
            """
        ))
    end

    for i in 1:D
        if length(axes[i].x) != size(data, i)
            throw(DimensionMismatch(
                """
                The `x` field in $(typeof(axes[i])) has length 
                $(length(axes[i].x)), but the `data` array has length 
                $(size(data, i)) at the corresponding dimension $i.
                These must match.
                """
            ))
        end
    end

    return NMRData{D, typeof(axes), typeof(data)}(axes, data)
end

export inv_out_1D
"""
    inv_out_1D(seq, x, y, xfit, yfit, X, f, r, SNR, α, selections, filter, title)

Output of the invert function for 1D pulse sequences.
A structure containing the following fields:
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PFG)
- `x`, the x values of the measurement (e.g time for relaxation or b-factor for diffusion).
- `y`, the y values of the measurement.
- `xfit`, the x values of the fitted data.
- `yfit`, the y values of the fitted data.
- `X`, the x values of the inversion results.
- `f`, the inversion results.
- `r`, the residuals.
- `SNR`, the signal-to-noise ratio.
- `α`, the regularization parameter.
- `selections`, a vector of tuples whose elements indicate the first and last index of the selected peaks.
- `filter`, a vector representing the mask used to filter and scale the data.
- `title`, a title describing the data.

"""
mutable struct inv_out_1D
    seq::Type{<:pulse_sequence1D}
    x::Vector
    y::Vector
    xfit::Vector
    yfit::Vector
    X::Vector
    f::Vector
    r::Vector
    SNR::Real
    alpha::Real
    selections::Vector{Tuple}
    filter::Vector
    title::String
end

export inv_out_2D
"""
    inv_out_2D(seq, X_direct, X_indirect, F, r, SNR, α, filter, selections)

Output of the invert function for 2D pulse sequences.
A structure containing the following fields:
- `seq` is the 2D pulse sequence (e.g. IRCPMG)
- `X_direct`, the x values of the direct dimension.
- `X_indirect`, the x values of the indirect dimension.
- `F`, the inversion results as a matrix.
- `r`, the residuals.
- `SNR`, the signal-to-noise ratio.
- `alpha`, the regularization parameter.
- `filter`, apply a mask to filter out artefacts when plotting.
- `selections`, the selection masks (e.g. when you want to highlight some peaks in a T₁-T₂ map).
- `title`, a title describing the data.

"""
mutable struct inv_out_2D
    seq::Type{<:pulse_sequence2D}
    X_direct::Vector
    X_indirect::Vector
    F::Matrix
    r::Vector
    SNR::Real
    alpha::Real
    filter::Matrix
    selections::Vector{Vector{Vector}}
    title::String
end


export expfit_struct
"""
Output of the expfit function.
Structure containing information about multiexponential fits.

The fields are as follows:
- `seq`: The pulse sequence
- `x` : The x acquisition values (e.g. time for relaxation or b-factor for diffusion).
- `y` : The y acquisition values.
- `xfit` : An array with fitted x values (for plotting purposes).
- `yfit` : An array with fitted y values (for plotting purposes).
- `u` : The fitted parameters for the `mexp` function.
- `u0` : The initial parameters for the `mexp` function.
- `r` : The residuals.
- `eq` : The equation of the fitted function.
- `eqn` : The equation of the initial function.
- `title` : A title describing the data.

"""
mutable struct expfit_struct
    seq::Type{<:NMRInversions.pulse_sequence1D}
    x::Vector
    y::Vector
    xfit::Vector
    yfit::Vector
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
