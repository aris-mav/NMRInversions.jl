## The following are custom types, defined mainly for multiple dispatch purposes

export pulse_sequence1D, pulse_sequence2D
export IR, SR, CPMG, PFG, IRCPMG, PFGCPMG, CPMGCPMG
export alpha_optimizer, regularization_solver

abstract type pulse_sequence1D end
abstract type pulse_sequence2D end

"Supported solvers for regularization"
abstract type regularization_solver end 

"Supported methods to determine the α parameter in tikhonov regularization"
abstract type alpha_optimizer end

"""
Inversion recovery pulse sequence for 1D relaxation experiments.
It can be used wherever the `seq` argument is required. 
"""
struct IR <: pulse_sequence1D end

"""
Saturation recovery pulse sequence for 1D relaxation experiments.
It can be used wherever the `seq` argument is required. 
"""
struct SR <: pulse_sequence1D end

"""
CPMG pulse sequence for 1D relaxation experiments.
It can be used wherever the `seq` argument is required. 
"""
struct CPMG <: pulse_sequence1D end

"""
Pulsed field gradient pulse sequence for 1D diffusion experiments.
It can be used wherever the `seq` argument is required. 
"""
struct PFG <: pulse_sequence1D end

"""
Inversion recovery - CPMG pulse sequence for 2D relaxation experiments (T1-T2).
The direct dimension is the T2, and the indirect dimension is the T1 acquisition times.
It can be used wherever the `seq` argument is required.
"""
struct IRCPMG <: pulse_sequence2D end

"""
PFG - CPMG pulse sequence for 2D D-T2 experiments.
It can be used wherever the `seq` argument is required.
"""
struct PFGCPMG <: pulse_sequence2D end

"""
CPMG - CPMG pulse sequence for 2D T2-T2 experiments.
It can be used wherever the `seq` argument is required.
"""
struct CPMGCPMG <: pulse_sequence2D end

export regularization_solver


export input1D
"""
    input1D(seq, x, y)
A structure containing the following fields:
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PFG)
- `x`, the x values of the measurement (e.g time for relaxation or b-factor for diffusion).
- `y`, the y values of the measurement. 

It can be used as an input for the invert and expfit functions.
"""
struct input1D
    seq::Type{<:pulse_sequence1D}
    x::AbstractVector{<:Real}
    y::AbstractVector
end

export input2D
"""
    input2D(seq, x, y)

A structure containing the following fields:
- `seq` is the 2D pulse sequence (e.g. IRCPMG)
- `x_direct` is the direct dimension acquisition parameter (e.g. the times when you aquire CPMG echoes).
- `x_indirect` is the indirect dimension acquisition parameter (e.g. all the delay times τ in your IR sequence).
- `data` is the 2D data matrix.

It can be used as an input for the invert function.
"""
struct input2D
    seq::Type{<:pulse_sequence2D}
    x_direct::AbstractVector{<:Real}
    x_indirect::AbstractVector{<:Real}
    data::AbstractMatrix
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
