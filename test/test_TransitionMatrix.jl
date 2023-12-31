using MarkovChainHammer, Test, LinearAlgebra
import MarkovChainHammer.TransitionMatrix: count_operator, generator, holding_times, perron_frobenius, sparse_count_operator
import MarkovChainHammer.TransitionMatrix: symmetric_generator, symmetric_perron_frobenius

@testset "Column Sum Consistency" begin
    timeseries = [1, 1, 1, 2, 2, 3, 3, 3, 2, 1]
    ones_vector = ones(3)
    zeros_vector = zeros(3)

    count_operator_computed = count_operator(timeseries)
    generator_computed = generator(timeseries)
    perron_frobenius_computed = perron_frobenius(timeseries)

    @test all([3, 3, 3] .== sum(count_operator_computed, dims=1)[:])
    @test all(zeros_vector .== sum(generator_computed, dims=1)[:])
    @test all(ones_vector .== sum(perron_frobenius_computed, dims=1)[:])
end

@testset "Hand Constructed: Default Case" begin
    # default case. See all states in timeseries
    timeseries = [1, 1, 1, 2, 2, 3, 3, 3, 2, 1]
    count_operator_exact = [2.0 1.0 0.0; 1.0 1.0 1.0; 0.0 1.0 2.0]
    generator_exact = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]
    perron_frobenius_exact = [2/3 1/3 0.0; 1/3 1/3 1/3; 0.0 1/3 2/3]
    perron_frobenius_exact_2step = [1/3 0 1/3; 2/3 0 1/3; 0 1 1/3]

    count_operator_computed = count_operator(timeseries)
    generator_computed = generator(timeseries)
    perron_frobenius_computed = perron_frobenius(timeseries)
    perron_frobenius_computed_2step = perron_frobenius(timeseries; step = 2)

    @test all(count_operator_exact .== count_operator_computed)
    @test all(generator_exact .== generator_computed)
    @test all(perron_frobenius_exact .== perron_frobenius_computed)
    @test all(perron_frobenius_exact_2step .== perron_frobenius_computed_2step)
end

@testset "Symmetric Hand Constructed: Default Case" begin
    # default case. See all states in timeseries
    timeseries = [1, 1, 1, 2, 2, 3, 3, 3, 2, 1]
    generator_exact = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]
    perron_frobenius_exact = [2/3 1/3 0.0; 1/3 1/3 1/3; 0.0 1/3 2/3]

    # trivial symmetry
    symmetries = [i -> i]

    generator_computed = symmetric_generator(timeseries, symmetries)
    perron_frobenius_computed = symmetric_perron_frobenius(timeseries, symmetries)
    @test all(generator_exact .== generator_computed)
    @test all(perron_frobenius_exact .== perron_frobenius_computed)

    # every state is equally likely symmetry
    function symmetry1(markov_index)
        if markov_index == 1
            return 2
        elseif markov_index == 2
            return 3
        elseif markov_index == 3
            return 1
        else
            return 0
        end
    end

    function symmetry2(markov_index)
        if markov_index == 1
            return 3
        elseif markov_index == 2
            return 1
        elseif markov_index == 3
            return 2
        else
            return 0
        end
    end

    symmetries = [symmetry1, symmetry2]

    perron_frobenius_exact = [5/9 2/9 2/9; 2/9 5/9 2/9; 2/9 2/9 5/9]

    generator_computed = symmetric_generator(timeseries, symmetries)
    perron_frobenius_computed = symmetric_perron_frobenius(timeseries, symmetries)

    @test all(generator_computed' .== generator_computed)
    @test all([generator_computed[i,i] .== generator_computed[1] for i in 1:3])

    @test all(perron_frobenius_computed' .== perron_frobenius_computed)
    @test all(perron_frobenius_exact .≈ perron_frobenius_computed)
    
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

@testset "Hand Constructed: Edge Case 2" begin
    # edge case. New state at the last step
    timeseries = [1, 1, 1, 2, 2, 4, 4, 4, 2, 1, 3]
    count_operator_exact = [2.0 1.0 0.0 0.0; 1.0 1.0 0.0 1.0; 1.0 0.0 0.0 0.0; 0.0 1.0 0.0 2.0]
    generator_exact = [-1/2 1/3 0.0 0.0; 1/4 -2/3 0.0 1/3; 1/4 0.0 0.0 0.0; 0.0 1/3 0.0 -1/3]
    perron_frobenius_exact = [1/2 1/3 0.0 0.0; 1/4 1/3 0.0 1/3; 1/4 0.0 1.0 0.0; 0.0 1/3 0.0 2/3]

    count_operator_computed = count_operator(timeseries)
    generator_computed = generator(timeseries)
    perron_frobenius_computed = perron_frobenius(timeseries)

    @test all(count_operator_exact .== count_operator_computed)
    @test all(generator_exact .== generator_computed)
    @test all(perron_frobenius_exact .== perron_frobenius_computed)
end

@testset "Hand Constructed: Edge Case 3 and Consistency" begin
    # edge case. New state at the last step
    timeseries = [1, 1, 1, 2, 2, 3, 3, 3, 2, 1, 4]
    permutation_matrix = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0] # reuse previously constructed matrix
    count_operator_exact = [2.0 1.0 0.0 0.0; 1.0 1.0 0.0 1.0; 1.0 0.0 0.0 0.0; 0.0 1.0 0.0 2.0]
    count_operator_exact = permutation_matrix * count_operator_exact * permutation_matrix'
    generator_exact = [-1/2 1/3 0.0 0.0; 1/4 -2/3 0.0 1/3; 1/4 0.0 0.0 0.0; 0.0 1/3 0.0 -1/3]
    generator_exact = permutation_matrix * generator_exact * permutation_matrix'
    perron_frobenius_exact = [1/2 1/3 0.0 0.0; 1/4 1/3 0.0 1/3; 1/4 0.0 1.0 0.0; 0.0 1/3 0.0 2/3]
    perron_frobenius_exact = permutation_matrix * perron_frobenius_exact * permutation_matrix'

    count_operator_computed = count_operator(timeseries)
    generator_computed = generator(timeseries)
    perron_frobenius_computed = perron_frobenius(timeseries)

    @test all(count_operator_exact .== count_operator_computed)
    @test all(generator_exact .== generator_computed)
    @test all(perron_frobenius_exact .== perron_frobenius_computed)
end

@testset "Transition Rate Matrix Consistency: Edge Case" begin
    timeseries = [1, 1, 1, 2, 2, 4, 4, 4, 2, 1]
    generator_computed1 = generator(timeseries)
    generator_computed2 = generator(timeseries, 4)
    generator_computed3 = generator(timeseries, 4, dt=1)
    generator_computed4 = generator(timeseries, 4, dt=0.5)

    @test all(generator_computed1 .== generator_computed2)
    @test all(generator_computed2 .== generator_computed3)
    @test all(generator_computed1 .== (0.5 * generator_computed4))
end

@testset "Hand Constructed: Default Case Sparse Check" begin
    # default case. See all states in timeseries
    timeseries = [1, 1, 1, 2, 2, 3, 3, 3, 2, 1]
    count_operator_exact = [2.0 1.0 0.0; 1.0 1.0 1.0; 0.0 1.0 2.0]
    generator_exact = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]
    perron_frobenius_exact = [2/3 1/3 0.0; 1/3 1/3 1/3; 0.0 1/3 2/3]
    perron_frobenius_exact_2step = [1/3 0 1/3; 2/3 0 1/3; 0 1 1/3]

    count_operator_computed = sparse_count_operator(timeseries)
    generator_computed = sparse_generator(timeseries)
    perron_frobenius_computed = sparse_perron_frobenius(timeseries)
    perron_frobenius_computed_2step = sparse_perron_frobenius(timeseries; step = 2)

    @test all(count_operator_exact .== count_operator_computed)
    @test all(generator_exact .== generator_computed)
    @test all(perron_frobenius_exact .== perron_frobenius_computed)
    @test all(perron_frobenius_exact_2step .== perron_frobenius_computed_2step)
end
