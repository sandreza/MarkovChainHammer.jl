using MarkovChainHammer, Test, Revise, Random, LinearAlgebra
import MarkovChainHammer.Utils: histogram
import MarkovChainHammer.Utils: autocovariance

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

