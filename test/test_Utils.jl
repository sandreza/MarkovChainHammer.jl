using MarkovChainHammer, Test, Revise, Random, LinearAlgebra
import MarkovChainHammer.Utils: histogram

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