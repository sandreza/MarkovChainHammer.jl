using MarkovChainHammer, Test, Random, LinearAlgebra
using MarkovChainHammer.Utils
import MarkovChainHammer.Utils: histogram
import MarkovChainHammer.Utils: autocovariance

@testset "utilities: scaled_entropy" begin
    for N in [10, 100, 123]
        p = ones(N) ./ N
        @test scaled_entropy(p) ≈ 1.0
    end
    p = [1.0, 0.0, 0.0]
    @test scaled_entropy(p) ≈ 0.0
    p = [0.5, 0.5, 0.0]
    @test scaled_entropy(p) ≈ -log(0.5) / log(3)
end

@testset "utilities: steady state" begin
    Q = [-1.0 1.0; 1.0 -1.0]
    p = steady_state(Q)
    @test all(p .≈ [0.5, 0.5])
    p = steady_state(exp(Q))
    @test all(p .≈ [0.5, 0.5])
end

@testset "utilities: koopman_modes" begin
    Q = [-1.0 2.0; 1.0 -2.0]
    p = steady_state(Q)
    km = koopman_modes(Q)
    @test abs(km[:, 1]' * p) < 100 * eps(1.0)
end

@testset "utilities: decorrelation times" begin
    Q = [-1.0 2.0; 1.0 -2.0]
    Λ = eigvals(Q)
    D = decorrelation_times(Q)
    @test D[end] == Inf
    @test D[1] == 1 / Λ[1]
end

@testset "Histogram Test: uniform weight" begin
    timeserieslist = [1 2 2 3 3 3 4 4 4 4 5 5 5 5 5]
    bins = 5
    xs, ys = histogram(timeserieslist, bins=bins)
    down, up = extrema(timeserieslist)
    vertices = range(down, up, length=bins + 1)
    bincenters = collect(vertices[2:end] + vertices[1:(end-1)]) .* 0.5
    @test length(xs) == bins
    @test length(ys) == bins
    @test all(xs .== bincenters)
    @test all(ys .== ys[1] .* sort(union(timeserieslist[:])))
end

@testset "Histogram Test: nonuniform weight" begin
    timeserieslist = [1 2 3 4 5 5]
    bins = 5
    p = [0.1, 0.2, 0.3, 0.4, 0.25, 0.25]
    xs, ys = histogram(timeserieslist, normalization=p, bins=bins)
    down, up = extrema(timeserieslist)
    vertices = range(down, up, length=bins + 1)
    bincenters = collect(vertices[2:end] + vertices[1:(end-1)]) .* 0.5
    @test length(xs) == bins
    @test length(ys) == bins
    @test all(xs .== bincenters)
    @test all(ys .≈ 0.1 .* collect(1:5))
end

@testset "Autocovariance Test" begin
    timeseries = [1 2 3 4 5 6 7 8 9 1]
    # timeseries check
    ac = autocovariance(timeseries)
    @test length(ac) == 10
    ac2 = autocovariance(timeseries; progress = true)
    @test length(ac2) == 10
    @test all((ac - ac2) .== 0)
    ac = autocovariance(timeseries; timesteps = 5)
    @test length(ac) == 5
    Q = generator(timeseries)
    # generator check
    g⃗ = ones(9)
    ac = autocovariance(g⃗, Q, 1:10)
    @test all(ac .< eps(100.0))
    Ps = [exp(Q * i) for i in 1:10]
    ac = autocovariance(g⃗, Ps, 10)
    @test all(ac .< eps(100.0))
end

@testset "Special Matrices" begin 
    Q = ornstein_uhlenbeck_generator(3)
    Qexact = [-1.0   0.5   0.0; 1.0  -1.0   1.0; 0.0   0.5  -1.0]
    @test all((Q .- Qexact) .≈ 0.0)
    A = central_advection_periodic(4)
    A_exact = -[  0.0   0.5   0.0  -0.5; -0.5   0.0   0.5   0.0; 0.0  -0.5   0.0   0.5; 0.5   0.0  -0.5   0.0]
    @test all((A .- A_exact) .≈ 0.0)
    Δ = discrete_laplacian_periodic(3)
    Δ_exact = [-2.0 1.0 1.0; 1.0 -2.0 1.0; 1.0 1.0 -2.0]
    @test all((Δ .- Δ_exact) .≈ 0.0)
end

@testset "Decomposition Test" begin
    for N in [3, 5, 7, 10]
        Q = ornstein_uhlenbeck_generator(N)
        E, R = exit_rate(Q)
        @test all(Q .≈ E * R)
        @test all([R[i,i] .≈ (N-1)/2 for i in 1:N])
        Qᴿ, Qⱽ = decomposition(Q)
        @test all(Qᴿ .≈ Q)
        @test norm(Qⱽ) ≤ eps(100.0 * N^2)
        @test all(Q .≈ Qᴿ + Qⱽ)
    end

    for N in [3, 5, 7, 10]
        advection = central_advection_periodic(N)
        diffusion = discrete_laplacian_periodic(N)
        Q = advection + diffusion
        E, R = exit_rate(Q)
        @test all(Q .≈ E * R)
        Qᴿ, Qⱽ = decomposition(Q)
        @test all(Qᴿ .≈ diffusion)
        @test all(norm(Qⱽ - advection) ≤ eps(100.0))
    end
end

