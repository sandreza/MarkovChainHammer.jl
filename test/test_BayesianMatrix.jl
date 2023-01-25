using MarkovChainHammer, LinearAlgebra, Test
using MarkovChainHammer.BayesianMatrix
using MarkovChainHammer.TransitionMatrix: generator
using Distributions, Random, Statistics

@testset "Bayesian Matrix: Basic Functionality" begin
    markov_chain = [1 1 1 1 2 2 1 1 1 3 1 1]

    Q = BayesianGenerator(markov_chain)
    @test rand(Q) isa Matrix
    @test length(Q.prior.rates) == 3
    @test length(Q.posterior.rates) == 3
    @test length(Q.prior.exit_probabilities) == 3
    @test length(Q.posterior.exit_probabilities) == 3

    Q̅ = mean(Q)
    Q_d = generator(markov_chain)

    @test all(abs.(Q̅ .- Q_d) .< eps(100.0))
    @test var(Q.posterior) isa Matrix
    @test all(size(Q) .== (3, 3))

    @test all(abs.(sqrt.(var(Q.posterior)) - std(Q.posterior) ) .< eps(100.0) )
    @test all(abs.(sqrt.(var(Q)) - std(Q) ) .< eps(100.0) )
    @test all(std(Q) - std(Q.posterior) .< eps(100.0))


    Q_r = rand(Q)
    Qs = rand(Q, 10)
    @test all(sum(Qs[3], dims=1) .< eps(100.0))


    ΛV = eigen(Q)
    @test ΛV isa Eigen

    Λ = eigvals(Q)
    @test Λ isa Vector

    ΛVs = eigen_distribution(Q; samples=10)
    @test length(ΛVs[1].values) == 3
    @test length(ΛVs) == 10
    @test ΛVs[1] isa Eigen

    Λs = eigvals_distribution(Q; samples=10)
    @test Λs isa Vector
end

@testset "Bayesian Matrix: Random Matrix Properties" begin
    markov_chain = [1 1 1 1 2 2 1 1 1 3 1 1]
    Q = BayesianGenerator(markov_chain)
    for N in [10, 100, 1000]
        Qs = rand(Q, N)
        μQ_empirical = mean(Qs)
        σ²Q_empirical = var(Qs)
        μQ = mean(Q)
        σ²Q = var(Q)
        @test all(abs.(μQ_empirical .- μQ) .< 10 / sqrt(N))
        @test all(abs.(σ²Q_empirical .- σ²Q) .< 10 / sqrt(N))
    end
end


@testset "Bayesian Matrix: Prior to Posterior" begin
    # test confidence with respect to same data
    markov_chain = [1 1 1 1 2 2 1 1 1 3 1 1]
    prior = GeneratorParameterDistributions(3)
    Q1 = BayesianGenerator(prior)
    Q2 = BayesianGenerator(markov_chain, prior)
    Q3 = BayesianGenerator(markov_chain, Q2.posterior)
    @test all(var(Q1) .> var(Q2))
    @test all(var(Q2) .> var(Q3))
    # test ordering
    markov_chain1 = [1 1 1 1 2 2 1 1 1 3 1 1]
    markov_chain2 = [1 2 1 2 1 3 1 3 2 2 1]

    prior = GeneratorParameterDistributions(3)
    Q1 = BayesianGenerator(markov_chain1, prior)
    Q2 = BayesianGenerator(markov_chain2, prior)
    Q21 = BayesianGenerator(markov_chain2, Q1.posterior)
    Q12 = BayesianGenerator(markov_chain1, Q2.posterior)

    @test sum([sum(abs.((params.(Q21.posterior.exit_probabilities)[i].-params.(Q12.posterior.exit_probabilities)[i])[1])) for i in 1:3]) .< eps(10.0)
    @test sum([sum(abs.((params.(Q21.posterior.rates)[i].-params.(Q12.posterior.rates)[i])[1])) for i in 1:3]) < eps(10.0)
end

