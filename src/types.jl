## The following are custom types for multiple dispatch purposes

# Pulse sequences
abstract type pulse_sequence1D end
abstract type pulse_sequence2D end

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


struct PFGCPMG <: pulse_sequence2D end
export pulse_sequence1D, pulse_sequence2D, IR, SR, CPMG, PFG, IRCPMG, PFGCPMG

# Supported solvers 
abstract type regularization_solver end 

"""
    brd
Solver for tikhonov (L2) regularization, following [this paper](https://doi.org/10.1109/78.995059)
from Venkataramanan et al.
Very fast, but only uses the identity as tiknohonov matrix.
No options required, it just works.
It can be used as a "solver" for the invert function.
"""
struct brd <: regularization_solver end


struct ripqp <: regularization_solver end

"""
    pdhgm(σ, τ)
Primal dual hybrid gradient method for L1 regularization, 
following [this paper](https://doi.org/10.1016/j.jmr.2017.05.010)
from Reci et al. / Journal of Magnetic Resonance 281 (2017) 188–198
It can be used as a "solver" for the invert function.

The particular choice of σ and τ is heuristic. 
A smaller σ will increase the stability while reducing the convergence speed
of the algorithm. A good compromise between the two was found when σ = 0.1 and τ = 10. 
The best values of σ and τ will depend slightly on the scaling of the signal. 
Therefore, it is best to normalize the NMR signal to a maximum of 1,
a technique which was followed in the cited study.
"""
struct pdhgm <: regularization_solver
    σ::Real
    τ::Real
end


"""
    optim_nnls(order)
Simple non-negative least squares method for tikhonov (L2) regularization, 
implemented using OptimizationOptimJl.
All around effective, but can be slow for large problems, such as 2D inversions.
It can be used as a "solver" for invert function.
Order is an integer that determines the tikhonov matrix 
(for more info look Hansen's 2010 book on inverse problems). 
Order `n` means that the penalty term will be the n'th derivative
of the results. 
"""
struct optim_nnls <: regularization_solver 
    order::Int
end

"""
    jump_nnls(order, solver)
Jump non-negative least squares method for tikhonov (L2) regularization, 
implemented using the JuMP extension.
All around effective, but can be slow for large problems, such as 2D inversions, 
unless a powerful solver like gurobi is used.
The solver argument is parsed directly into JuMP.
It can be used as a "solver" for invert function.
Order determines the tikhonov matrix. If 0 is chosen, the identity matrix is used.
"""
struct jump_nnls <: regularization_solver 
    order::Int
    solver::Symbol
end

export regularization_solver, brd, ripqp, pdhgm, optim_nnls, jump_nnls

# Supported methods to determine regularization α parameter
abstract type smoothing_optimizer end

"""
    gcv
Generalized cross validation for finding the optimal regularization parameter α.
"""
struct gcv <: smoothing_optimizer end


"""
    lcurve(a,b,n)
L-curve method for finding the optimal regularization parameter α.
It will test a set of logarithmically-spaced values, 
starting from `a`, ending in `b`, and consisting of `n` values.
The optimal value will be chosen based on the maximum curvature of the L curve,
as described in Hansen 2010, "Discrete Inverse Problems".

This is usually less reliable than the `gcv` method, 
but it's nice to have multiple options.

Note that the first time you use this in every session, 
it will take about a minute to compile the code.
This might be optimized in the future, if there's demand for it.

"""
struct lcurve <: smoothing_optimizer 
    lowest_value::Real
    highest_value::Real
    number_of_steps::Int
end

export smoothing_optimizer, gcv, lcurve


export svd_kernel_struct
"""
    svd_kernel_struct(K,g,U,S,V)
A structure containing the following fields:
- `K`, the kernel matrix.
- `G`, the data vector.
- `U`, the left singular values matrix.
- `S`, the singular values vector.
- `V`, the right singular values matrix.

To access the fields of a structure, we use the dot notation, 
e.g. if the structure is called `a` and we want to access the kernel contained in there,
we type `a.K`
"""
struct svd_kernel_struct
    K::AbstractMatrix
    g::AbstractVector
    U::AbstractMatrix
    S::AbstractVector
    V::AbstractMatrix
end


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
    inv_out_1D(seq, x, y, xfit, yfit, X, f, r, SNR, α, selections, title)

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
- `selections`, a vector of tuples whose elements are indicate the first and last index of the selected peaks. 
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
    title::String
end



export inv_out_2D
"""
    inv_out_2D(seq, X_dir, X_indir, F, r, SNR, α, filter, selections)

Output of the invert function for 2D pulse sequences.
A structure containing the following fields:
- `seq` is the 2D pulse sequence (e.g. IRCPMG)
- `X_dir`, the x values of the direct dimension.
- `X_indir`, the x values of the indirect dimension.
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
    X_dir::Vector
    X_indir::Vector
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
