using Test
using NMRInversions
using SparseArrays
using LinearAlgebra
using Random
using GLMakie

Random.seed!(1)

function test_artificial_data(seq::Type{<:pulse_sequence1D}, SNR = 100 ; kwargs...)

    x = exp10.(range(log10(1e-4), log10(5), 32)) # acquisition range
    X = exp10.(range(-5, 1, 128)) # T range
    K = create_kernel(seq, x, X)
    f_custom = [0.5exp.(-(x)^2 / 3) + exp.(-(x - 1.3)^2 / 0.5) for x in range(-5, 5, length(X))]
    g = K * f_custom
    y = g +  (maximum(g)/SNR) .* randn(length(x))
    results = invert( seq, x, y, lims=(-5,1,128);kwargs...) 

    score = LinearAlgebra.norm(f_custom - (results.f ./ (maximum(results.f)) .* maximum(f_custom)))
    display("Sequence: $seq, SNR: $SNR, Score: $score")
    
    # The condition below is COMPLETELY arbitrary
    # just happens to match what we expect as sensible results
    return score < 10 * SNR^(-0.35)

    #=f = Figure()=#
    #=ax = Axis(f[:,:])=#
    #=plot!(ax, f_custom )=#
    #=plot!(ax, results.f ./ (maximum(results.f)) .* maximum(f_custom))=#
    #=return f=#
end

function test_artificial_data_2D() #throws errors, needs fixing

    x_direct = exp10.(range(log10(1e-4), log10(5), 1024)) # acquisition range
    x_indirect = exp10.(range(log10(1e-4), log10(5), 32)) # acquisition range

    X_direct = exp10.(range(-5, 1, 64)) # T range
    X_indirect = exp10.(range(-5, 1, 64)) # T range

    # Create a rotated gaussian distribution (which would look like real data)
    θ = 135
    σ₁ = 1.3
    σ₂ = 0.4
    x₀ = 0
    y₀ = 1.3
    a = ((cosd(θ)^2) / (2 * σ₁^2)) + ((sind(θ)^2) / (2 * σ₂^2))
    b = -((sind(2 * θ)) / (4 * σ₁^2)) + ((sind(2 * θ)) / (4 * σ₂^2))
    c = ((sind(θ)^2) / (2 * σ₁^2)) + ((cosd(θ)^2) / (2 * σ₂^2))
    F_original = ([exp.(-(a * (x - x₀)^2 + 2 * b * (x - x₀) * (y - y₀) + c * (y - y₀)^2)) for x in range(-5, 5, length(X_direct)), y in range(-5, 5, length(X_indirect))])

    K1 = create_kernel(CPMG, x_direct, X_direct)
    K2 = create_kernel(IR, x_indirect, X_indirect)

    data = K1 * F_original * K2'
    data = complex.(data, 0.001 .* maximum(real(data)) .* randn(size(data)))

    results = invert(IRCPMG, x_direct, x_indirect, data, alpha=0.01, lims1=(-5, 1, 64), lims2=(-5, 1, 64), normalize=false)

    score = LinearAlgebra.norm(results.F - F_original)

    return  score < 0.5

end


function test_phase_correction(plots=false)

    # Create real and imaginary parts
    Re_original = exp.(-range(1, 20, 1000)) + randn(1000) .* 0.01
    Im_original = randn(1000) .* 0.01

    # Get them out of phase
    ϕd = rand() * 2π
    Re_shifted, Im_shifted = NMRInversions.phase_shift(Re_original, Im_original, ϕd)

    # Correct the phase
    Rₙ, Iₙ, ϕc = NMRInversions.autophase(Re_shifted, Im_shifted, 1)

    ## Plots for sanity check (using Plots.jl)
    if plots == true
        p1 = plot([Re_original, Im_original], label=["Original real" "Original Imaginary"])
        p2 = plot([Re_shifted, Im_shifted], label=["Dephased real" "Dephased Imaginary"])
        p3 = plot([Rₙ, Iₙ], label=["Corrected real" "Corrected Imaginary"])
        ϕ_range = range(0, 2π, 20000)
        Re1_vs_φ = Re_shifted[1] .* cos.(ϕ_range) - Im_shifted[1] .* sin.(ϕ_range)
        Im_sum_vs_φ = [im_cost([ϕ], (Re_shifted, Im_shifted)) for ϕ in ϕ_range]
        p4 = plot(ϕ_range, Re1_vs_φ, xlabel="ϕ", label="Re[1]")
        p4 = plot!(ϕ_range, (Im_sum_vs_φ ./ maximum(Im_sum_vs_φ)) .* maximum(Re1_vs_φ), xlabel="ϕ", label="sum(im.^2)", legend=:topleft)
        p4 = vline!(p4, [ϕc], label="corrected phase")
        display(plot(p1, p2, p3, p4))
    end

    display("The correction error is $(round(2π - (ϕd + ϕc), sigdigits=2)) radians")

    return abs(2π - (ϕd + ϕc)) < 0.05
end


function test_expfit()

    x = [range(0.001, 3, 32)...]
    u = [3, 0.2 ,4 ,0.05]
    y = mexp(CPMG, u, x) + 0.01 .* randn(length(x))
    data = input1D(CPMG, x, y)
    results = expfit(2,data,normalize=false)

    return sum((sort(u) .- sort(results.u)) .^ 2) < 0.5
end


@testset "NMRInversions.jl" begin

    display("Testing inversion of artificial data.")
    for seq in [IR, CPMG]
        for snr in exp10.(2:1:4) 
            @test test_artificial_data(seq,snr,alpha=gcv(),silent=true)
        end
    end

    @test test_artificial_data(IR, 500, alpha=gcv(1))
    @test test_artificial_data(IR, 500, alpha=lcurve(1e-5,10))

    display("Testing expfit.")
    @test test_expfit()

    display("Testing phase correction accuracy.")
    @test test_phase_correction()

    display("Testing plots - inverting geospec data.")
    @test isa( 
        plot(invert(import_geospec("../example_data/geospec_data/bunter_IR.txt"),silent=true)), 
        Figure
    )
    @test isa( 
        plot(invert(
            import_geospec("../example_data/geospec_data/bunter_IRCPMG.txt"),
            alpha = 0.5, silent = true)
             ), 
        Figure
    )

end
