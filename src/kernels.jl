export create_kernel

export svd_kernel_struct
"""
    svd_kernel_struct(K,g,U,S,V)
A structure containing the following fields:
- `K`, the kernel matrix.
- `g`, the data vector.
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
    _kernel_eq(axis<:DataAxis, X; x0, y, n)

Internal functions passed to the one defined below.

Arguments:

- `axis` is the `DataAxis` for the experiment.
- `X` is the array for T1, T2, D etc.. (can also be a DataAxis)
- `x0` is the offset in `x` (can be the first point `x[1]`).
- `y` is the recorded data, used to normalise the kernel.
- `n` should be `1` for exponential and `2` for gaussian decays.
"""
_kernel_eq(axis::IR, X; x0, y, n) =
    @. y[end] - (y[end] + abs(y[1])) * exp(-((axis.x - x0) / X)^n)

_kernel_eq(axis::SR, X; x0, y, n) =
    @. y[end] * (1 - exp(-((axis.x - x0) / X))^n)

_kernel_eq(axis::CPMG, X; x0, y, n) =
    @. y[1] * exp(-((axis.x - x0) / X)^n)

_kernel_eq(axis::PFG, X; x0, y, n) =
    @. y[1] * exp(-((axis.x - x0) * X)^n)


"""
    create_kernel(x, X)

Arguments:
- `x` is the `DataAxis` (x-axis) for the experiment (time, b-factor etc.)
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
    x::DataAxis,
    X::AbstractVector{<:Real};
    y::AbstractVector=ones(1),
    gaussian::Bool=false,
    x0::Real=0
)
    return _kernel_eq(
        x, X';
        y=real.(y),
        n=gaussian ? 2 : 1,
        x0=(isa(x, CPMG) ? x0 : x.x[1]),
    )
end


"""
1D kernel with SVD compression.
"""
function create_kernel(
    input::ExperimentData{1},
    X::NTuple{1,AbstractVector},
)

    W_½ = sqrt(input.W[1])

    axis = input.axes[1]
    g = W_½ * input.data

    usv = svd(W_½ * create_kernel(axis, X[1], y=g))

    idx = isnan(input.SNR) ? (:) : findall(i -> i .> (1 / input.SNR), usv.S)

    U = usv.U[:, idx]
    S = usv.S[idx]
    V = usv.V[:, idx]

    K_new = Diagonal(S) * V'
    g_new = U' * real.(g)

    return svd_kernel_struct(K_new, g_new, U, S, V)
end


"""
2D kernel with SVD compression (following Mitchell 2012 paper).
"""
function create_kernel(
    input::ExperimentData{2},
    X::NTuple{2,AbstractVector},
)

    axes = input.axes
    data = input.data

    W1_½ = sqrt(input.W[1])
    W2_½ = sqrt(input.W[2])

    G = W1_½ * real.(data) * W2_½

    # Generate Kernels
    K_dir = create_kernel(axes[1], X[1], y=vec(data[:, end]))
    K_indir = create_kernel(axes[2], X[2], y=vec(data[1, :]))

    ## Perform SVD truncation
    usv_dir = svd(W1_½ * K_dir) #paper (13)
    usv_indir = svd(W2_½ * K_indir) #paper (14)

    # finding which singular components are contributed from K1 and K2
    # using outer product instead of Kronecker, to make indices more obvious
    S21 = usv_indir.S * usv_dir.S'

    idx = isnan(input.SNR) ? (:) : findall(i -> i .> (1 / input.SNR), S21)

    s̃ = S21[idx]
    ñ = length(s̃)

    si = (first.(Tuple.(idx)))  # direct dimension singular vector indices
    sj = (last.(Tuple.(idx)))   # indirect dimension singular vector indices

    g̃ = diag(usv_indir.U[:, si]' * G' * usv_dir.U[:, sj])
    Ũ₀ = Array{Float64}(undef, 0, 0) # no such thing as U in this case, it is absorbed by g 

    Ṽ₀ = zeros(size(usv_indir.V, 1) * size(usv_dir.V, 1), ñ)
    for k in 1:ñ
        # Get the specific vectors for this significant component
        v_indir = usv_indir.V[:, si[k]]
        v_dir = usv_dir.V[:, sj[k]]

        # The Kronecker product combines them into a single basis vector
        # This represents the joint "space" of both dimensions
        Ṽ₀[:, k] = kron(v_indir, v_dir)
    end

    K̃₀ = Diagonal(s̃) * Ṽ₀'

    return svd_kernel_struct(K̃₀, g̃, Ũ₀, s̃, Ṽ₀)
end
