using MarkovChainHammer, Test, Revise
import MarkovChainHammer.TransitionMatrix: count_operator, generator, holding_times, perron_frobenius

@testset "Hand Constructed: Default Case" begin
    # default case. See all states in timeseries
    timeseries = [1, 1, 1, 2, 2, 3, 3, 3, 2, 1]
    count_operator_exact = [2.0 1.0 0.0; 1.0 1.0 1.0; 0.0 1.0 2.0]
    generator_exact = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]
    perron_frobenius_exact = [2/3 1/3 0.0; 1/3 1/3 1/3; 0.0 1/3 2/3]

    count_operator_computed = count_operator(timeseries)
    generator_computed = generator(timeseries)
    perron_frobenius_computed = perron_frobenius(timeseries)

    @test all(count_operator_exact .== count_operator_computed)
    @test all(generator_exact .== generator_computed)
    @test all(perron_frobenius_exact .== perron_frobenius_computed)
end

@testset "Hand Constructed: Edge Case" begin
    # edge case. Only see 3 out of 4 states
    timeseries = [1, 1, 1, 2, 2, 4, 4, 4, 2, 1]
    count_operator_exact = [2.0 1.0 0.0 0.0; 1.0 1.0 0.0 1.0; 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 2.0]
    generator_exact = [-1/2 1/3 0.0 0.0; 1/2 -2/3 0.0 1/3; 0.0 0.0 0.0 0.0; 0.0 1/3 0.0 -1/3]
    perron_frobenius_exact = [2/3 1/3 0.0 0.0; 1/3 1/3 0.0 1/3; 0.0 0.0 1.0 0.0; 0.0 1/3 0.0 2/3]

    count_operator_computed = count_operator(timeseries)
    generator_computed = generator(timeseries)
    perron_frobenius_computed = perron_frobenius(timeseries)

    @test all(count_operator_exact .== count_operator_computed)
    @test all(generator_exact .== generator_computed)
    @test all(perron_frobenius_exact .== perron_frobenius_computed)
end

@testset "Transition Rate Matrix Consistenty: Edge Case" begin
    timeseries = [1, 1, 1, 2, 2, 4, 4, 4, 2, 1]
    generator_computed1 = generator(timeseries)
    generator_computed2 = generator(timeseries, 4)
    generator_computed3 = generator(timeseries, 4, dt = 1)
    generator_computed4 = generator(timeseries, 4, dt = 0.5)

    @test all(generator_computed1 .== generator_computed2)
    @test all(generator_computed2 .== generator_computed3)
    @test all(generator_computed1 .==  (0.5 * generator_computed4))
end
