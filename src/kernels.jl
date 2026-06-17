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

- `axis` is the `DataAxis` for the experiment.
- `X` is the array for T1, T2, D etc..

The ones below are set as keyword arguments so that they 
are not broadcasted in `create_kernel`.

- `x0` is the offset in `x` (can be the first point `x[1]`).
- `y` is the recorded data, used to normalise the kernel.
- `n` should be `1` for exponential and `2` for gaussian decays.
"""
kernel_eq(axis::Type{IR}, X; x0, y, n) = 
    y[end] - (y[end] + abs(y[1])) * exp(-((axis.x-x0) / X)^n)

kernel_eq(axis::Type{SR}, X; x0, y, n) = 
    y[end] * (1 - exp(-((axis.x-x0) / X))^n)

kernel_eq(axis::Type{CPMG}, X; x0, y, n) = 
    y[1] * exp(-((axis.x-x0) / X)^n)

kernel_eq(axis::Type{PFG}, X; x0, y, n) = 
    y[1] * exp(-((axis.x-x0) * X)^n)


"""
# Create a kernel for the inversion of 1D data.
    create_kernel(axis, X)

Arguments:
- `axis` is the `DataAxis` for the experiment.
- `X` is the vector for the output's x-axis (T1, T2, D etc.)

Keyword (optional) arguments:

- `y`, if you provide the measurement data `y`, the kernel 
will be normalised to match the highest value in `y`. Defaults 
to an array containing `1.0`.
- `gaussian` determines whether the kernel is exponential or 
gaussian. Defaults to `false` (exponential by default).
- `x0` is the offset in `x` (can be the first point `x[1]`).

The output is a matrix, `K`.

"""
function create_kernel(
    axis::DataAxis, X::Vector; 
    y::Vector=ones(1), 
    gaussian::Bool = false, 
    x0::Real=0
)
    return kernel_eq.(
        axis.x, X'; 
        y= real.(y), 
        n= gaussian ? 2 : 1,
        x0= (isa(axis, CPMG) ? x0 : axis.x[1]), 
    )
end


"""
    create_kernel(seq, x, X, g)
If data vector of real values is provided, SVD is performed on the kernel, and the output is a "svd_kernel_struct" instead.

If data vector is complex, the SNR is calculated and the SVD is automatically truncated accordingly,
to remove the "noisy" singular values.
"""
function create_kernel(axis::DataAxis, X::Vector, g::Vector{<:Real})

    usv = svd(create_kernel(axis , X, y=g))
    K_new = Diagonal(usv.S) * usv.V'
    g_new = usv.U' * g

    return svd_kernel_struct(K_new, g_new, usv.U, usv.S, usv.V)

end


function create_kernel(axis::DataAxis, X::Vector, g::Vector{<:Complex})

    usv = svd(create_kernel(axis , X, y=g))

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


## Multidimensional cases

"Passing the 1D tuple to the one of the 1D functions."
function create_kernel(
    axes::NTuple{1, DataAxis},
    X::NTuple{1, AbstractVector},
    Data::AbstractMatrix{<:Complex};
    kwargs...
)
    create_kernel(axes[1], X[1], Data; kwargs...)
end


"""
# Generating a kernel for a 2D inversion
    create_kernel(seq, x_direct, x_indirect, X_direct, X_indirect, Data)

- `axis` is a tuple of the two `DataAxis` objects.
- `X` is a tuple of the two vectors containing the output values (T1,T2,D etc..).
- `Data` is the 2D data matrix of complex data.
The output is an `svd_kernel_struct`

"""
function create_kernel(
    axes::NTuple{2, DataAxis},
    X::NTuple{2, AbstractVector},
    Data::AbstractMatrix{<:Complex}
)

    G = real.(Data)
    SNR = calc_snr(Data)

    # Generate Kernels
    K_dir = create_kernel(axes[1], X[1], y = vec(Data[:,end]))
    K_indir = create_kernel(axes[2], X[2], y = vec(Data[1,:]))

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
