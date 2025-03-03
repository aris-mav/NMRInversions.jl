
"""
# Create a kernel for the inversion of 1D data.
    create_kernel(seq, x, X)
- `seq` is the pulse sequence (e.g. IR, CPMG, PFG)
- `x` is the experiment x axis (time or b factor etc.)
- `X` is the range for the output x axis (T1, T2, D etc.)

The output is a matrix, `K`.

"""
function create_kernel(seq::Type{<:pulse_sequence1D}, x::Vector, X::Vector)

    if seq == IR
        kernel_eq = (t, T) -> 1 - 2 * exp(-t / T)
    elseif seq == SR
        kernel_eq = (t, T) -> 1 - exp(-t / T)
    elseif seq == CPMG
        kernel_eq = (t, T) -> exp(-t / T)
    elseif seq == PFG
        kernel_eq = (b, D) -> exp(-b * D)
    end

    return kernel_eq.(x, X')
end

"""
    create_kernel(seq, x, X, g)
If data vector of real values is provided, SVD is performed on the kernel, and the output is a "svd_kernel_struct" instead.

If data vector is complex, the SNR is calculated and the SVD is automatically truncated accordingly,
to remove the "noisy" singular values.
"""
function create_kernel(seq::Type{<:pulse_sequence1D}, x::Vector, X::Vector, g::Vector{<:Real})

    usv = svd(create_kernel(seq, x , X))
    K_new = Diagonal(usv.S) * usv.V'
    g_new = usv.U' * g

    return svd_kernel_struct(K_new, g_new, usv.U, usv.S, usv.V)

end

function create_kernel(seq::Type{<:pulse_sequence1D}, x::Vector, X::Vector, g::Vector{<:Complex})

    usv = svd(create_kernel(seq, x , X))

    SNR = calc_snr(g)
    indices = findall(i -> i .> (1 / SNR), usv.S) # find elements in S12 above the noise threshold

    display("SVD truncated to $(length(indices)) singular values out of $(length(usv.S))")

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

    if SNR < 1000
        @warn("The SNR is $(round(SNR, digits=1)), which is below the recommended value of 1000. 
            Consider running experiment with more scans.")
    end

    # Generate Kernels
    if seq == IRCPMG
        K_dir = create_kernel(CPMG, x_direct, X_direct)
        K_indir = create_kernel(IR, x_indirect, X_indirect)
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
    Ṽ₀ = repeat(usv_dir.V[:, sj], size(usv_indir.V, 1), 1) .* reshape(repeat(usv_indir.V[:, si]', size(usv_dir.V, 1), 1), ñ, size(usv_indir.V,1)*size(usv_dir.V,1) )'
    K̃₀ = Diagonal(s̃) * Ṽ₀'

    return svd_kernel_struct(K̃₀, g̃, Ũ₀, s̃, Ṽ₀)
end

