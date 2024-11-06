"""
Create a finite difference differentiation matrix of order `order` for a matrix of size `m x m`.
"""
function Γ(m::Int, order::Int)
    # Eilers, P. H. C. (2003). Analytical Chemistry, 75(14), 3631–3636. (Supporting Information)
    # Taken from RegularizationTools.jl
    if order == 0
        return Array{Float64}(LinearAlgebra.I, (m, m))
    end
    return diff(Γ(m, order - 1), dims=1)
end

"""
    calc_snr(data)
Calculate the Signal-to-Noise Ratio (SNR) from complex data,
where the real part is (mostly) signal and the imaginary part is (mostly) noise.
The STD of the latter half of the imaginary signal is used for the calculation 
(former half might contain signal residues as well). 
"""
function calc_snr(data::AbstractMatrix{<:Complex}) # For matrix data

    real_data = real.(data)
    imag_data = imag.(data)

    noise = imag_data[(floor(Int64, (size(imag_data, 1) / 2))):end, :]
    σ_n = sqrt(sum((noise .- sum(noise) / length(noise)) .^ 2) / (length(noise) - 1))
    SNR = maximum(abs.(real_data)) / σ_n

    return SNR
end
function calc_snr(data::AbstractVector{<:Complex}) # For vector data

    real_data = real.(data)
    imag_data = imag.(data)

    noise = imag_data[(floor(Int64, (size(imag_data, 1) / 2))):end]
    σ_n = sqrt(sum((noise .- sum(noise) / length(noise)) .^ 2) / (length(noise) - 1))
    SNR = maximum(abs.(real_data)) / σ_n

    return SNR
end


"""
Return a vector of matrices, containing the F for each selection polygon.
"""
function selections(res::inv_out_2D)

    dir = res.X_dir
    indir = res.X_indir
    F = res.F

    x = collect(1:length(indir))
    y = collect(1:length(dir))

    z = res.F' .* res.filter'

    points = [[i, j] for i in x, j in y] #important, used by inpolygon later
    mask = zeros(size(points))

    polygons = res.selections 
    selections = [zeros(size(F)) for _ in 1:length(polygons)]

    for (i, polygon) in enumerate(polygons)
        mask .= [PolygonOps.inpolygon(p, polygon; in=1, on=1, out=0) for p in points]
        selections[i] .= mask .* z 
    end

    return selections

end

function calc_weighted_averages(F::Matrix, dir::Vector, indir::Vector)

    T1 = vec(sum(F, dims=2)) ⋅ indir / sum(F)
    T2 = vec(sum(F, dims=1)) ⋅ dir / sum(F)

    return [T1, T2]
end


function entropy(u::AbstractArray, p::Tuple)

    phc0 = u[1]
    phc1 = u[2]
    R⁰ = p[1]
    I⁰ = p[2]
    γ = p[3]

    n = length(R⁰)
    ϕ = phc0 .+ phc1 .* collect(range(1, n) ./ n)

    R, _ = phase_shift(R⁰, I⁰, ϕ)

    h = abs.(diff(R)) ./ sum(abs.(diff(R)))

    return -sum(h .* log.(h)) + γ * sum((R .^ 2)[findall(x -> x < 0, R)])

end


function minimize_entropy(Re, Im, γ)

    optf = Optimization.OptimizationFunction(entropy)
    prob = Optimization.OptimizationProblem(optf, [0.1, 0.1], (Re, Im, γ))# lb=[0.0001, 0.0001], ub=[2π, 2π])

    uopt = OptimizationOptimJL.solve(prob, OptimizationOptimJL.SimulatedAnnealing())

    n = length(Re)
    ϕ = uopt[1] .+ uopt[2] .* collect(range(1, n) ./ n)

    Rₙ, Iₙ = phase_shift(Re, Im, ϕ)

    p = plot(a[:, 1], Rₙ, label="Real")
    p = plot!(a[:, 1], Iₙ, label="Imaginary")
    display(p)

    return Rₙ, Iₙ

end


"""
Cost function - Sum of squares of the imaginary Part of NMR data.
"""
function im_cost(u, p)

    Re = p[1]
    Im = p[2]
    ϕ = u[1]

    _, Iₙ = phase_shift(Re, Im, ϕ)
    return sum(Iₙ .^ 2)
end

"""
Shift the phase of a complex signal by φ radians.
Re is the real part, and Im is the imaginary.
"""
function phase_shift(Re, Im, ϕ)

    Re_new = Re .* cos(ϕ) - Im .* sin(ϕ)
    Im_new = Im .* cos(ϕ) + Re .* sin(ϕ)

    return Re_new, Im_new
end


"""
# Auto phase correcton for input data.
This function corrects the phase of the inputs so that the real part is mostly pure data and the imaginary part is mostly noise.
The output is a structure of the same type as the input, but the data within will be phase-corrected.
"""
function autophase(input::input1D)
    re = real.(input.y)
    im = imag.(input.y)
    seq = input.seq

    if seq in [IR]
        re, im, ϕ = autophase(re, im, -1)
        return input1D(seq, input.x, complex.(re, im))

    elseif seq in [CPMG, PFG]
        re, im, ϕ = autophase(re, im, 1)
        return input1D(seq, input.x, complex.(re, im))
    else
        error("This function is not yet implemented for this pulse sequence. Please submit an issue on GitHub.")
    end

end

function autophase(input::input2D)
    re = real.(input.data)
    im = imag.(input.data)
    seq = input.seq

    if seq in [IRCPMG]
        re, im, ϕ = autophase(re, im, -1)
        return input2D(seq, input.x_direct, input.x_indirect, complex.(re, im))
    else
        error("This function is not yet implemented for this pulse sequence. Please submit an issue on GitHub.")
    end

end

"""
## Full syntax:.
    autophase(re, im, startingpoint)
- `re` is the real part of the signal.
- `im` is the imaginary part of the signal.
- `startingpoint` is a value between -1 and 1. It basically determines what should the 1st point of the real part be (e.g. -1 for IR, 1 for CPMG).
"""
function autophase(re, im, startingpoint::Real) # startingpoint is a value between -1 and 1

    ϕ_range = range(0, 2π, 500)
    Re1_vs_φ = re[1] .* cos.(ϕ_range) - im[1] .* sin.(ϕ_range)

    ϕ₀ = ϕ_range[argmin(abs.(Re1_vs_φ .- startingpoint * maximum(Re1_vs_φ)))]

    optf = Optimization.OptimizationFunction(im_cost, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(optf, [ϕ₀], (re, im), lb=[ϕ₀ - 1], ub=[ϕ₀ + 1], x_tol=1e-5)
    ϕ = OptimizationOptimJL.solve(prob, OptimizationOptimJL.BFGS())[1]

    Rₙ, Iₙ = phase_shift(re, im, ϕ)

    return Rₙ, Iₙ, ϕ
end


## Data compression (NOT WORKING YET)

function compress_data(t_direct::AbstractVector, G::AbstractMatrix, bins::Int=64)
    # Compress direct dimension to the length of bins

    # Use logarithmically spaced windows for a window average
    windows = zeros(bins)

    x = 1
    while windows[1] == 0
        windows = (exp10.(range(log10(x), log10(10), bins))) # make log array
        windows = (windows / sum(windows)) * size(G)[1] # normalize it so that elements add up to correct size
        windows = floor.(Int, LowerTriangular(ones(Int, bins, bins)) * windows) # make valid indices
        x += 1
        #sanity check, this sum below should be almost equal to the uncompressed array size
        #sum(windows[2:end]-windows[1:end-1]) 
    end

    W0 = Diagonal(inv(LowerTriangular(ones(Int, bins, bins))) * windows)

    # Window average matrix, A
    A = zeros(bins, size(G)[1])
    for i in 1:bins
        if i == 1
            A[i, 1:windows[1]] .= 1 / diag(W0)[i]
        else
            A[i, (windows[i-1]+1):windows[i]] .= 1 / diag(W0)[i]
        end
    end

    t_direct = A * t_direct # Replace old direct time array with compressed one
    G = A * G # Replace old G with compressed one

    # sanity check plot:
    usv1 = svd(sqrt(W0) * K1) #paper (13)
    # surface(G, camera=(110, 25), xscale=:log10)

end


export weighted_averages
"""
    weighted_averages(r::inv_out_1D)
Return a vector with the weighted averages for the selections in `r.selections`.
"""
function weighted_averages(r::inv_out_1D)

    wa = Vector(undef, length(r.selections))
    for (i,s) in enumerate(r.selections)
        wa[i] = r.f[s[1]:s[2]]' * r.X[s[1]:s[2]] / sum(r.f[s[1]:s[2]]) 
    end
    return wa
end
