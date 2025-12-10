"""
# Inversion for 1D pulse sequences:
    invert(seq, x, y ; lims, alpha, solver, normalize)

This function will build a kernel and use it to perform an inversion using the algorithm of your choice.
The output is an `inv_out_1D` structure.

 Necessary (positional) arguments:
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PFG)
- `x` is the experiment x axis (time or b factor etc.)
- `y` is the experiment y axis (intensity of the NMR signal)

 Optional (keyword) arguments:
- `lims=(a,b,c)` will set the "limits" of the output X, 
so that it starts from 10^a, ends in 10^b and consists of c 
logarithmically spaced values.
Alternatively, a vector of values can be used directly, if more freedom is needed 
(e.g. `lims=exp10.(range(-4, 1, 128))`).
By default, 128 points are used, logarithmically spaced between the first value 
in x divided by 7, and the final value in x multiplied by 7.

- `alpha` determines the smoothing term. Use a real number for a fixed alpha, 
or look at the "Types and Structures - Finding Optimal Alpha" section of the 
documentation for other options.  
No selection will lead to automatically determining alpha through the 
default method, which is `gcv()`.
- `solver` is the algorithm used to do the inversion math. Default is `brd`.
- `normalize` will normalize `y` to 1 at the max value of `y`. Default is `true`.  

"""
function invert(seq::Type{<:pulse_sequence1D}, x::AbstractArray, y::Vector;
                lims::Union{Tuple{Real, Real, Int}, AbstractVector, Type{<:pulse_sequence1D}}=seq, 
                alpha::Union{Real, alpha_optimizer}=gcv(), 
                solver::Union{regularization_solver, Type{<:regularization_solver}}=brd(),
                normalize::Bool = true,
                silent::Bool = true
                )

    if normalize
        y = y ./ y[argmax(real(y))]
    end

    if isa(lims, Tuple)
        X = exp10.(range(lims...))
    elseif isa(lims, AbstractVector)
        X = lims
    elseif isa(lims, Type{<:pulse_sequence1D})
        if lims == PFG
            X = collect(NMRInversions.logrange(minimum(x)*7 , maximum(x)*50 , 128)) * 1e-9
        else
            X = collect(NMRInversions.logrange(minimum(x)/7 , maximum(x)*7 , 128)) 
        end
    end


    # Change scale to match bfactor, which is s/m²e-9, undo below to go back to SI
    if seq == PFG
        X .= X .* 1e9
    end
 
    ker_struct = create_kernel(seq, x, X, y)
    α = 0.0 #placeholder, will be replaced below 

    if isa(alpha, Real)

        α = alpha
        f, _ = solve_regularization(ker_struct.K, ker_struct.g, α, solver)

    else

        f, _, α = find_alpha(ker_struct, solver, alpha ; silent = silent)

    end

    x_fit = exp10.(range(log10(x[1]), log10(x[end]), 512))
    y_fit = create_kernel(seq, x_fit, X, y=y) * f

    isreal(y) ? SNR = NaN : SNR = calc_snr(y)
    
    r = create_kernel(seq, x, X, y=y) * f - y

    if seq == PFG
        X .= X ./ 1e9
    end

    return inv_out_1D(seq, x, y, x_fit, y_fit, X, f, r, SNR, α, [], ones(length(X)),"")

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
    invert(seq, x_direct, x_indirect, Data; lims1, lims2, alpha, solver, normalize)

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

In both cases above, you can use a tuple specifying the limits of the range, 
or a vector of values, same as the `lims` argument in the 1D inversion function.
By default, 64 points are used, logarithmically spaced between the first value 
in x divided by 7, and the final value in x multiplied by 7 (for both direct and indirect x).

- `alpha` determines the smoothing term. Use a real number for a fixed alpha.  No selection will lead to automatically determining alpha through the default method, which is `gcv`.
- `solver` is the algorithm used to do the inversion math. Default is `brd`.
- `normalize` will normalize `Data` so that its maximum value is 1. Default is `true`.

"""
function invert(
    seq::Type{<:pulse_sequence2D}, x_direct::AbstractVector, x_indirect::AbstractVector, Data::AbstractMatrix;
    lims1::Union{Tuple{Real, Real, Int}, AbstractVector, Type{<:pulse_sequence2D}}=seq, 
    lims2::Union{Tuple{Real, Real, Int}, AbstractVector, Type{<:pulse_sequence2D}}=seq,
    alpha::Union{Real, alpha_optimizer} = gcv(), 
    solver::Union{regularization_solver, Type{<:regularization_solver}}=brd(),
    normalize::Bool=true,
    silent::Bool = false
)

    if normalize
        Data = Data ./ Data[argmax(real(Data))]
    end

    if isa(lims1, Type{<:pulse_sequence2D})
        X_direct = collect(NMRInversions.logrange(minimum(x_direct)/7 , maximum(x_direct)*7 , 64)) 
    elseif isa(lims1, Tuple)
        X_direct = exp10.(range(lims1...))
    elseif isa(lims1, AbstractVector)
        X_direct = lims1
    end

    if isa(lims2, Type{<:pulse_sequence2D})
        if lims2 == IRCPMG
            X_indirect = collect(NMRInversions.logrange(minimum(x_indirect)/7, maximum(x_indirect)*7, 64)) 
        elseif lims2 == PFGCPMG
            X_indirect = 1e-18 .* collect(NMRInversions.logrange(minimum(x_indirect)*7, maximum(x_indirect)*50, 64)) 
        end
    elseif isa(lims2, Tuple)
        X_indirect = exp10.(range(lims2...))
    elseif isa(lims2, AbstractVector)
        X_indirect = lims2
    end


    ker_struct = create_kernel(seq, x_direct, x_indirect, X_direct, X_indirect, Data)

    α = 1.0 #placeholder, will be replaced below

    if isa(alpha, Real)
        α = alpha
        f, r = solve_regularization(ker_struct.K, ker_struct.g, α, solver)

    else
        f, r, α = find_alpha(ker_struct, solver, alpha, silent = silent)

    end

    F = reshape(f, length(X_direct), length(X_indirect))

    results = inv_out_2D(seq, X_direct, X_indirect, F, r, calc_snr(Data), α, ones(size(F)), [],"")

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
