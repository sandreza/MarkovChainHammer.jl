using MarkovChainHammer, Test, Revise, Random, LinearAlgebra
import MarkovChainHammer.TransitionMatrix: count_operator, generator, holding_times, perron_frobenius
import MarkovChainHammer.Trajectory: generate
Random.seed!(12345)

@testset "Convergence Test: Generator" begin
    errorlist = []
    errorlist2 = []
    Q_exact = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]

    for dt in [0.3, 1.0, 3.0]
        P_exact = exp(dt * Q_exact)
        for n in [10^i for i in 1:6]
            markov_chain = generate(Q_exact, n; dt=dt)
            Q = generator(markov_chain, 3; dt=dt)
            P = perron_frobenius(markov_chain, 3)
            push!(errorlist, norm(P - P_exact) / norm(P_exact))
            @test(errorlist[end] < 6 / sqrt(n))
            P2 = exp(Q * dt)
            push!(errorlist2, norm(P2 - P_exact) / norm(P_exact))
            push!(plist, P)
            self_consistent_error = 2 * norm(P2 - P) / norm(P)
            @test(errorlist2[end] < maximum([6 / sqrt(n), self_consistent_error]))
        end
    end
end

@testset "Convergence Test: Perron Frobenius" begin
    errorlist = []
    errorlist2 = []
    Q_exact = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]

    for dt in [0.3, 1.0, 3.0]
        P_exact = exp(dt * Q_exact)
        for n in [10^i for i in 1:6]
            markov_chain = generate(P_exact, n; dt=dt) # different from previous
            Q = generator(markov_chain, 3; dt=dt)
            P = perron_frobenius(markov_chain, 3)
            push!(errorlist, norm(P - P_exact) / norm(P_exact))
            @test(errorlist[end] < 6 / sqrt(n))
            P2 = exp(Q * dt)
            push!(errorlist2, norm(P2 - P_exact) / norm(P_exact))
            push!(plist, P)
            self_consistent_error = 2 * norm(P2 - P) / norm(P)
            @test(errorlist2[end] < maximum([6 / sqrt(n), self_consistent_error]))
        end
    end
end

@testset "Convergence Test: Heuristic" begin
    Q_exact = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]
    default_size = length(generate(Q))
    @test(default_size <= length(generate(0.1 * Q)))
    @test(default_size >= length(generate(10 * Q)))
    @test(default_size >= length(generate(100 * Q)))
    @test(2 == length(generate(0 * Q)))

end