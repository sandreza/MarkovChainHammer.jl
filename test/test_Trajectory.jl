using MarkovChainHammer, Test, Revise, Random, LinearAlgebra
import MarkovChainHammer.TransitionMatrix: count_operator, generator, holding_times, perron_frobenius
import MarkovChainHammer.Trajectory: generate
Random.seed!(12345)

@testset "Convergence Test: Generator" begin
    Q_exact = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]
    for dt in [0.3, 1.0, 3.0]
        P_exact = exp(dt * Q_exact)
        for n in [10^i for i in 1:6]
            markov_chain = generate(Q_exact, n; dt=dt)
            Q = generator(markov_chain, 3; dt=dt)
            P = perron_frobenius(markov_chain, 3)
            perron_frobenius_error = norm(P - P_exact) / norm(P_exact)
            @test(perron_frobenius_error < 6 / sqrt(n))
            P2 = exp(Q * dt)
            perron_frobenius_error2 = norm(P2 - P_exact) / norm(P_exact)
            self_consistent_error = 2 * norm(P2 - P) / norm(P)
            @test(perron_frobenius_error2 < maximum([6 / sqrt(n), self_consistent_error]))
        end
    end
end

@testset "Convergence Test: Perron Frobenius" begin
    Q_exact = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]
    for dt in [0.3, 1.0, 3.0]
        P_exact = exp(dt * Q_exact)
        for n in [10^i for i in 1:6]
            markov_chain = generate(P_exact, n; dt=dt) # different from previous
            Q = generator(markov_chain, 3; dt=dt)
            P = perron_frobenius(markov_chain, 3)
            perron_frobenius_error = norm(P - P_exact) / norm(P_exact)
            @test(perron_frobenius_error < 6 / sqrt(n))
            P2 = exp(Q * dt)
            perron_frobenius_error2 = norm(P2 - P_exact) / norm(P_exact)
            self_consistent_error = 2 * norm(P2 - P) / norm(P)
            @test(perron_frobenius_error2 < maximum([6 / sqrt(n), self_consistent_error]))
        end
    end
end

@testset "Step Heuristic Test" begin
    Q = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]
    default_size = length(generate(Q))
    @test(default_size <= length(generate(0.1 * Q)))
    @test(default_size >= length(generate(10 * Q)))
    @test(default_size >= length(generate(100 * Q)))
    @test(2 == length(generate(0 * Q)))
end


@testset "Holding Times" begin
    markov_chain = Int.([ones(5); ones(4) * 2; ones(2)])
    markov_chain_holding_times = holding_times(markov_chain)
    @test all(markov_chain_holding_times[1] .== [5, 2])
    @test all(markov_chain_holding_times[2] .== 4)

    markov_chain_holding_times = holding_times(markov_chain; dt=3.0)

    @test all(markov_chain_holding_times[1] .== [15.0, 6.0])
    @test all(markov_chain_holding_times[2] .== 12.0)
end