"""
# Inversion for 1D pulse sequences:
    invert(seq, x, y ; lims, alpha, solver)

This function will build a kernel and use it to perform an inversion using the algorithm of your choice.
The output is an `inv_out_1D` structure.

 Necessary (positional) arguments:
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PGSE)
- `x` is the experiment x axis (time or b factor etc.)
- `y` is the experiment y axis (intensity of the NMR signal)

 Optional (keyword) arguments:
- `lims=(a,b,c)` will set the "limits" of the output X, 
so that it starts from 10^a, ends in 10^b and consists of c 
logarithmically spaced values 
(default is (-5, 1, 128) for relaxation and (-10, -7, 128) for diffusion). 
Alternatiively, a vector of values can be used directly, if more freedom is needed 
(e.g. `lims=exp10.(range(-5, 1, 128))`).
- `alpha` determines the smoothing term. Use a real number for a fixed alpha.  No selection will lead to automatically determining alpha through the defeault method, which is `gcv`.
- `solver` is the algorithm used to do the inversion math. Default is `brd`.

"""
function invert(seq::Type{<:pulse_sequence1D}, x::AbstractArray, y::Vector;
    lims::Union{Tuple{Real, Real, Int}, AbstractVector, Type{<:pulse_sequence1D}}=seq, 
    alpha::Union{Real, smoothing_optimizer, Type{<:smoothing_optimizer}}=gcv, 
    solver::Union{regularization_solver, Type{<:regularization_solver}}=brd)

    if isa(lims, Tuple)
        X = exp10.(range(lims...))
    elseif isa(lims, AbstractVector)
        X = lims
    elseif isa(lims, Type{<:pulse_sequence1D})
        if lims == PFG
            X = exp10.(range(-10, -7, 128))
        else
            X = exp10.(range(-5, 1, 128))
        end
    end

    α = 1 #placeholder, will be replaced below 

    if isa(alpha, Real)

        α = alpha
        ker_struct = create_kernel(seq, x, X, y)

        f, r = solve_regularization(ker_struct.K, ker_struct.g, α, solver)

    elseif alpha == gcv

        ker_struct = create_kernel(seq, x, X, y)
        f, r, α = solve_gcv(ker_struct, solver)

    else
        error("alpha must be a real number or gcv")

    end

    x_fit = collect(range(0, x[end] + 0.1 * x[end], 128))
    y_fit = create_kernel(seq, x_fit, X) * f

    isreal(y) ? SNR = NaN : SNR = calc_snr(y)

    weighted_average =  f'X / sum(f)

    return inv_out_1D(seq, x, y, x_fit, y_fit, X, f, r, SNR, α, weighted_average)

end

"""
    invert(data::input1D ; kwargs...)

Instead of the positional arguments `seq`, `x` and `y`,
you can use a single `input1D` structure, which contains the same information. 
Especially useful if you're using the output of one 
of the import functions (look documentation tutorial section).
"""
function invert(data::input1D; kwargs...)

    return invert(data.seq, data.x, data.y; kwargs...)

end



"""
# Inversion for 2D pulse sequences:
    invert(seq, x_direct, x_indirect, X_direct, X_indirect, Data)

This function will build a kernel and use it to perform an inversion using the algorithm of your choice.
The output is an `inv_out_2D` structure.

 Necessary (positional) arguments:
- `seq` is the 2D pulse sequence (e.g. IRCPMG)
- `x_direct` is the direct dimension acquisition parameter (e.g. the times when you aquire CPMG echoes).
- `x_indirect` is the indirect dimension acquisition parameter (e.g. all the delay times τ in your IR sequence).
- `Data` is the 2D data matrix of complex data.


 Optional (keyword) arguments:
- `lims1` determines the output "range" of the inversion in the direct dimension (e.g. T₂ times in IRCPMG)
- `lims2` determines the output "range" of the inversion in the indirect dimension (e.g. T₁ times in IRCPMG)
 In both cases above, you can use a tuple specifying the limits of the range, or a vector of values, same as the `lims` argument in the 1D inversion.

- `alpha` determines the smoothing term. Use a real number for a fixed alpha.  No selection will lead to automatically determining alpha through the defeault method, which is `gcv`.
- `solver` is the algorithm used to do the inversion math. Default is `brd`.

"""
function invert(
    seq::Type{<:pulse_sequence2D}, x_direct::AbstractVector, x_indirect::AbstractVector, Data::AbstractMatrix;
    lims1::Union{Tuple{Real, Real, Int}, AbstractVector}=(-5, 1, 100), 
    lims2::Union{Tuple{Real, Real, Int}, AbstractVector}=(-5, 1, 100),
    alpha::Union{Real, smoothing_optimizer, Type{<:smoothing_optimizer}}=gcv, 
    solver::Union{regularization_solver, Type{<:regularization_solver}}=brd)

    if isa(lims1, Tuple)
        X_direct = exp10.(range(lims1...))
    elseif isa(lims1, AbstractVector)
        X_direct = lims1
    end

    if isa(lims2, Tuple)
        X_indirect = exp10.(range(lims2...))
    elseif isa(lims2, AbstractVector)
        X_indirect = lims2
    end


    ker_struct = create_kernel(seq, x_direct, x_indirect, X_direct, X_indirect, Data)

    α = 1.0 #placeholder, will be replaced below

    if isa(alpha, Real)
        α = alpha
        f, r = solve_regularization(ker_struct.K, ker_struct.g, α, solver)

    elseif alpha == gcv
        f, r, α = solve_gcv(ker_struct, solver)

    else
        error("alpha must be a real number or gcv")

    end

    F = reshape(f, length(X_direct), length(X_indirect))

    results = inv_out_2D(seq, X_direct, X_indirect, F, r, calc_snr(Data), α, ones(size(F)), [])

    return results

end


"""
    invert(data::input2D ; kwargs...)

Instead of the positional arguments `seq`, `x_direct` ,
`x_indirect` and `Data`, you can use a single `input2D` 
structure, which contains the same information. 
Especially useful if you're using the output of one of the 
import functions (look documentation tutorial section).
"""
function invert(input::input2D; kwargs...)

    invert(input.seq, input.x_direct, input.x_indirect, input.data; kwargs...)

end
