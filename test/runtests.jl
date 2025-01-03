using Test
using NMRInversions
using SparseArrays
using LinearAlgebra


seq = IR
x = exp10.(range(log10(1e-4), log10(5), 32)) # acquisition range
X = exp10.(range(-5, 1, 128)) # T range
K = create_kernel(seq, x, X)
f_custom = [0.5exp.(-(x)^2 / 3) + exp.(-(x - 1.3)^2 / 0.5) for x in range(-5, 5, length(X))]
g = K * f_custom
y = g + 0.001 * maximum(g) .* randn(length(x))


function test_lcurve()

    x = collect(range(0.01, 1, 32))
    X = exp10.(range(-5, 1, 128)) # T range
    K = create_kernel(CPMG, x, X)
    f_custom = [0.5exp.(-(x)^2 / 3) + exp.(-(x - 1.3)^2 / 0.5) for x in range(-5, 5, length(X))]
    g = K * f_custom
    noise_level = 0.001 * maximum(g)
    y = g + noise_level .* randn(length(x))

    alphas = exp10.(range(log10(1e-5), log10(1e0), 128))
    curvatures = zeros(length(alphas))
    xis = zeros(length(alphas))
    rhos = zeros(length(alphas))

    for (i, α) in enumerate(alphas)
        println("α = ", α)

        f, r = NMRInversions.solve_regularization(K, y, α, brd)

        ξ = f'f
        ρ = r'r
        λ = √α

        A = sparse([K; √(α) * LinearAlgebra.I ])
        b = sparse([r; zeros(size(A, 1) - size(r, 1))])

        z = NMRInversions.solve_ls(A, b)

        ∂ξ∂λ = (4 / λ) * f'z

        ĉ = 2 * (ξ * ρ / ∂ξ∂λ) * (α * ∂ξ∂λ * ρ + 2 * ξ * λ * ρ + λ^4 * ξ * ∂ξ∂λ) / ((α * ξ^2 + ρ^2)^(3 / 2))

        xis[i] = ξ
        rhos[i] = ρ
        curvatures[i] = ĉ

    end

    α = alphas[argmin(curvatures)]

    f, r = NMRInversions.solve_regularization(K, y, α, brd)

    #=using Plots=#
    #=begin=#
    #=    p1 = plot(alphas, curvatures, xscale=:log10, xlabel="α", ylabel="curvature",label = "curvature vs. α");=#
    #=    p1 = vline!(p1, [α], label="α = $α");=#
    #=    p2 = plot(X, [f_custom, f], label=["original" "solution"], xscale=:log10);=#
    #=    p3 = plot(rhos, xis, xscale=:log10, yscale=:log10,label = "lcurve");=#
    #=    p3 = scatter!([rhos[argmin(curvatures)]], [xis[argmin(curvatures)]], label="α = $α");=#
    #=    p4 = scatter(x, y, label="data");=#
    #=    p4 = plot!(x, K * f, label="solution");=#
    #=    p = plot(p1, p2, p3, p4)=#
    #=end=#
end

function testT1T2()

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

    return LinearAlgebra.norm(results.F - F_original) < 0.5

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
    @test norm(invert(
        seq, x, y, alpha=gcv(), lims=(-5,1,128), normalize=false
    ).f - f_custom) < 1.0
    @test norm(invert(
        seq, x, y, alpha=gcv(1), lims=(-5,1,128), normalize=false
    ).f - f_custom) < 1.0
    @test norm(invert(
        seq, x, y, alpha=lcurve(1e-5, 10), lims=(-5,1,128), normalize=false
    ).f - f_custom) < 1.0
    @test testT1T2()
    @test test_expfit()
    @test test_phase_correction()

end
