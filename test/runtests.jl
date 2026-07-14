using Test
using NMRInversions
using SparseArrays
using LinearAlgebra
using Random
using GLMakie

Random.seed!(1)

function test_invert(seq::Type{<:DataAxis}, SNR=100; kwargs...)

    x = seq(
        exp10.(range(log10(1e-4), log10(5), 128))
    ) # acquisition range
    X = exp10.(range(-5, 1, 128)) # T range
    K = create_kernel(x, X)
    f_custom = [
        0.5exp.(-(x)^2 / 3) + exp.(-(x - 1.3)^2 / 0.5)
        for x in range(-5, 5, length(X))
    ]
    g = K * f_custom
    y = @. g + (maximum(g) / SNR) * randn()

    data = ExperimentData(seq(x), y)
    results = invert(data, axes=(X,); kwargs...)

    # The condition below is COMPLETELY arbitrary, but seems sensible
    score = LinearAlgebra.norm(
        f_custom - (results.data ./ (maximum(results.data)) .* maximum(f_custom))
    )
    threshold = 10 * SNR^(-0.35)

    display("Sequence: $seq, SNR: $SNR, Score: $score, Threshold: $threshold")
    return score < threshold

end

function T2T2_data()

    x_direct = range(1, 5000, 500) .* (250 * 1e-6)
    x_indirect = logrange(1, 5000, 32) .* (250 * 1e-6)
    X = exp10.(range(-5, 1, 128))
    f = [
        0.5exp.(-(x + 0.7)^2 / 0.2) + exp.(-(x - 1.3)^2 / 1.9)
        for x in range(-5, 5, length(X))
    ]
    A = create_kernel(CPMG, x_direct, X)
    B = create_kernel(CPMG, x_indirect, X)
    K = Diagonal(f)

    G = A * K * B'

    noise_real = (maximum(G) * 0.01) .* randn(size(G)...)
    noise_imag = (maximum(G) * 0.01) .* randn(size(G)...)

    data = ExperimentData(
        (CPMG(x_direct), CPMG(x_indirect)),
        complex.(G .+ noise_real, noise_imag)
    )

    return data
end

function test_expfit(seq::Type{<:DataAxis})

    x = seq(collect(NMRInversions.logrange(0.001, 3, 32)))
    u = [3, 0.2, 4, 0.05]
    y = mexp(u, x) + 0.01 .* randn(length(x))
    data = ExperimentData(x, y)

    results = expfit(data, 2)

    return sum((sort(u) .- sort(results.u)) .^ 2) < 0.5
end

@testset "Inversions on artificial data" begin
    for seq in [IR, CPMG, SR]
        for snr in exp10.(2:1:4)
            @test test_invert(seq, snr, alpha=GCV(), silent=true)
        end
    end
    @test test_invert(IR, 500, alpha=GCV(1), silent=true)
    @test test_invert(IR, 500, alpha=LCurve(1), silent=true)
    @test test_invert(IR, 500, alpha=GCV(1e-5, 10), silent=true)
    @test test_invert(IR, 500, alpha=LCurve(1e-5, 10), silent=true)
end

@testset "expfits" begin
    for seq in [IR, CPMG, SR]
        @test test_expfit(seq)
    end
end

@testset "GLMakie and import_geospec" begin
    @test isa(
        plot(invert(
            import_geospec("../example_data/geospec_txt/IR_bunter.txt"),
            silent=true)
        ),
        Figure
    )
    @test isa(
        plot(invert(
            import_geospec("../example_data/geospec_txt/IRCPMG_bunter.txt"),
            alpha=0.5, silent=true)
        ),
        Figure
    )
end
