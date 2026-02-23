export create_kernel

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

"""
    kernel_eq
# Internal function
Return the value of the kernel at a given element.

Arguments:

- `seq` is the pulse sequence (e.g. IR, CPMG, PFG)
- `x` is the experiment x axis (time or b factor etc.)
- `X` is the range for the output x axis (T1, T2, D etc.)

The ones below are set as keyword arguments so that they 
are not broadcasted in `create_kernel`.

- `x0` is the offset in `x` (normally the first point `x[1]`).
- `y` is the recorded data, used to normalise the kernel.
- `n` should be `1` for exponential and `2` for gaussian decays.
"""
kernel_eq(::Type{IR}, x, X; x0, y, n) = y[end] - (y[end] + abs(y[1])) * exp(-((x-x0) / X)^n)
kernel_eq(::Type{SR}, x, X; x0, y, n) = y[end] * (1 - exp(-((x-x0) / X))^n)
kernel_eq(::Type{CPMG}, x, X; x0, y, n) = y[1] * exp(-((x-x0) / X)^n)
kernel_eq(::Type{PFG}, x, X; x0, y, n) = y[1] * exp(-((x-x0) * X)^n)


"""
# Create a kernel for the inversion of 1D data.
    create_kernel(seq, x, X)

Arguments:
- `seq` is the pulse sequence (e.g. IR, CPMG, PFG)
- `x` is the experiment x axis (time or b factor etc.)
- `X` is the range for the output x axis (T1, T2, D etc.)

Keyword (optional) arguments:

- `y`, if you provide the measurement data `y`, the kernel 
will be normalised to match the highest value in `y`. Defaults 
to an array containing `1.0`.
- `gaussian` determines whether the kernel is exponential or 
gaussian. Defaults to `false` (exponential by default).

The output is a matrix, `K`.

"""
function create_kernel(seq::Type{<:pulse_sequence1D}, x::Vector, X::Vector; 
                       y::Vector=ones(1), gaussian = false)

    return kernel_eq.(seq, x, X'; x0= x[1], y= real.(y), n= gaussian ? 2 : 1)

end


"""
    create_kernel(seq, x, X, g)
If data vector of real values is provided, SVD is performed on the kernel, and the output is a "svd_kernel_struct" instead.

If data vector is complex, the SNR is calculated and the SVD is automatically truncated accordingly,
to remove the "noisy" singular values.
"""
function create_kernel(seq::Type{<:pulse_sequence1D}, x::Vector, X::Vector, g::Vector{<:Real})

    usv = svd(create_kernel(seq, x , X, y=g))
    K_new = Diagonal(usv.S) * usv.V'
    g_new = usv.U' * g

    return svd_kernel_struct(K_new, g_new, usv.U, usv.S, usv.V)

end

function create_kernel(seq::Type{<:pulse_sequence1D}, x::Vector, X::Vector, g::Vector{<:Complex})

    usv = svd(create_kernel(seq, x , X, y = g))

    SNR = calc_snr(g)
    indices = findall(i -> i .> (1 / SNR), usv.S) # find elements in S12 above the noise threshold

    #=display("SVD truncated to $(length(indices)) singular values out of $(length(usv.S))")=#

    U = usv.U[:, indices]
    S = usv.S[indices]
    V = usv.V[:, indices]

    K_new = Diagonal(S) * V'
    g_new = U' * real.(g)

    return svd_kernel_struct(K_new, g_new, U, S, V)

end

"""
# Generating a kernel for a 2D inversion
    create_kernel(seq, x_direct, x_indirect, X_direct, X_indirect, Data)

- `seq` is the 2D pulse sequence (e.g. IRCPMG)
- `x_direct` is the direct dimension acquisition parameter (e.g. the times when you aquire CPMG echoes).
- `x_indirect` is the indirect dimension acquisition parameter (e.g. all the delay times τ in your IR sequence).
- `X_direct` is output "range" of the inversion in the direct dimension (e.g. T₂ times in IRCPMG)
- `X_indirect` is output "range" of the inversion in the indirect dimension (e.g. T₁ times in IRCPMG)
- `Data` is the 2D data matrix of complex data.
The output is an `svd_kernel_struct`

"""
function create_kernel(seq::Type{<:pulse_sequence2D},
    x_direct::AbstractVector, x_indirect::AbstractVector,
    X_direct::AbstractVector, X_indirect::AbstractVector,
    Data::AbstractMatrix{<:Complex})

    G = real.(Data)
    SNR = calc_snr(Data)

    # Generate Kernels
    if seq == IRCPMG
        K_dir = create_kernel(CPMG, x_direct, X_direct, y = vec(Data[:,end]))
        K_indir = create_kernel(IR, x_indirect, X_indirect, y = vec(Data[1,:]))
    elseif seq == PFGCPMG
        K_dir = create_kernel(CPMG, x_direct, X_direct, y = vec(Data[:,1]))
        K_indir = create_kernel(PFG, x_indirect, X_indirect, y = vec(Data[1,:]))
    elseif seq == CPMGCPMG
        K_dir = create_kernel(CPMG, x_direct, X_direct, y = vec(Data[:,1]))
        K_indir = create_kernel(CPMG, x_indirect, X_indirect, y = vec(Data[1,:]))
    else
        error("2D inversion not yet implemented for $(seq)")
    end

    ## Perform SVD truncation
    usv_dir = svd(K_dir) #paper (13)
    usv_indir = svd(K_indir) #paper (14)

    # finding which singular components are contributed from K1 and K2
    S21 = usv_indir.S * usv_dir.S' #Using outer product instead of Kronecker, to make indices more obvious
    indices = findall(i -> i .> (1 / SNR), S21) # find elements in S12 above the noise threshold
    s̃ = S21[indices]
    ñ = length(s̃)

    si = (first.(Tuple.(indices)))  # direct dimension singular vector indices
    sj = (last.(Tuple.(indices)))   # indirect dimension singular vector indices

    g̃ = diag(usv_indir.U[:, si]' * G' * usv_dir.U[:, sj])
    Ũ₀ = Array{Float64}(undef, 0, 0) # no such thing as U in this case, it is absorbed by g 

    Ṽ₀ = zeros(size(usv_indir.V, 1) * size(usv_dir.V, 1), ñ)
    for k in 1:ñ
        # Get the specific vectors for this significant component
        v_indir = usv_indir.V[:, si[k]]
        v_dir   = usv_dir.V[:, sj[k]]

        # The Kronecker product combines them into a single basis vector
        # This represents the joint "space" of both dimensions
        Ṽ₀[:, k] = kron(v_indir, v_dir)
    end

    K̃₀ = Diagonal(s̃) * Ṽ₀'

    return svd_kernel_struct(K̃₀, g̃, Ũ₀, s̃, Ṽ₀)
end

