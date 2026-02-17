using LinearAlgebra

"""
    hankel(x::Vector, L::Int)

Constructs a Hankel matrix from vector `x` with `L` rows.
No zero padding.
"""
function hankel(x::Vector{T}, L::Int) where T
    N = length(x)
    K = N - L + 1
    return [x[i + j - 1] for i in 1:L, j in 1:K]
end
"""
    hankel(x::Vector, L::Int)

Constructs a Hankel matrix from vector `x` with zero padding
at the end (MATLAB style).
"""
function hankel(x::Vector{T}) where T
    N = length(x)
    
    return [ (i + j - 1) <= N ? x[i + j - 1] : zero(T) for i in 1:N, j in 1:N]
end

function make_data()

    timepoints = 200;
    T2 = 1e-3; 
    dt = 1e-4; 
    amp = 1;

    time = collect(0:timepoints-1)*dt
    data = amp.* exp.( - time ./ T2) 

    return input1D(CPMG, time, data) 
end

function solve_MPM(data::input1D)

    f = real.(data.y)
    n = length(data)
    dt = data.x[2] - data.x[1]

    ncomp = 1
    datamatrix = hankel(f)
    u, s_vec, v = svd(datamatrix)
    s = Diagonal(s_vec)

    decaymatrix_svd = zeros(size(datamatrix)) 
    for ii = 1:ncomp
        decaymatrix_svd = u[:, ii] * s[ii, ii] * v[:, ii]'
    end

    n = length(data)
    d1 = zeros(n-1, n-1)
    d2 = zeros(n-1, n-1)

    for ii = 1:(n - 1)
        d1[ii, :] = decaymatrix_svd[ii, 1:(n - 1)]
        d2[ii, :] = decaymatrix_svd[ii, 2:n]
    end

    eigs = eigvals(pinv(d1) * d2)

    eigvec_left = zeros(n-1, ncomp)
    eigvec_right = zeros(ncomp, n-1)
    T_mpm = zeros(ncomp)

    for jj = 1:ncomp
        for ii = 1:(n - 1)
            eigvec_left[ii, jj] = real(eigs[jj])^(ii - 1)
            eigvec_right[jj, ii] = real(eigs[jj])^(ii - 1)
        end
        T_mpm[jj] = -dt / real(log(Complex(real(eigs[jj]))))
    end

    amplitude_mpm = pinv(eigvec_left) * decaymatrix_svd[1:(n-1), 1:(n-1)] * pinv(eigvec_right)

    return T_mpm, amplitude_mpm
end
