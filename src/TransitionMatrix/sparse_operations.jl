using SparseArrays, ProgressBars
export sparse_perron_frobenius, sparse_generator

function sparse_count_operator(markov_chain::Vector{S}, number_of_states::S, step::Int) where S
    count_matrix = spzeros(typeof(markov_chain[1]), number_of_states, number_of_states)
    for i in ProgressBar(0:step-1)
        reduced_markov_chain = markov_chain[1+i:step:end]
        for j in 1:length(reduced_markov_chain)-1
            count_matrix[reduced_markov_chain[j+1], reduced_markov_chain[j]] += 1
        end
    end
    return count_matrix
end

function sparse_count_operator(markov_chain::Vector{S}, number_of_states::S) where S
    @info "computing sparse count matrix"
    count_matrix = spzeros(S, number_of_states, number_of_states)
    for i in ProgressBar(1:length(markov_chain)-1)
        count_matrix[markov_chain[i+1], markov_chain[i]] += 1
    end
    return count_matrix
end

sparse_count_operator(markov_chain) = sparse_count_operator(markov_chain, maximum(markov_chain))

"""
sparse_perron_frobenius(markov_chain; step = 1)

# Description 
    Calculate the perron-frobenius matrix from a markov chain in sparse format

# Arguments
- `markov_chain::AbstractVector`: A vector of integers representing the state of a markov chain at each time step.

# Keyword Arguments
- `step::Integer=1`: The step size of the constructed operator.

# Returns
- `perron_frobenius_matrix::Matrix`: The perron-frobenius matrix of the markov chain in sparse format
"""
function sparse_perron_frobenius(partitions::Vector{S}; step = 1) where S
    number_of_states = maximum(partitions)
    count_matrix = sparse_count_operator(partitions, number_of_states, step)
    number_of_states = maximum(partitions)
    perron_frobenius_matrix = Float64.(count_matrix) 
    normalization = sum(count_matrix, dims=1)
    for i in ProgressBar(eachindex(normalization))
        for j in perron_frobenius_matrix[:, i].nzind
            perron_frobenius_matrix[j, i] /= normalization[i]
        end
    end
    return perron_frobenius_matrix
end

"""
sparse_generator(markov_chain; dt = 1)

# Description 
    Calculate the generator matrix from a markov chain in sparse format

# Arguments
- `markov_chain::AbstractVector`: A vector of integers representing the state of a markov chain at each time step.

# Keyword Arguments
- `dt::Real=1`: The time step of the constructed operator.

# Returns
- `generator_matrix::Matrix`: The generator matrix of the markov chain in sparse format
"""
function sparse_generator(partitions::Vector{S}; dt = 1) where S
    number_of_states = maximum(partitions)
    count_matrix = sparse_count_operator(partitions, number_of_states)
    generator_matrix = Float64.(count_matrix) 
    for i in 1:number_of_states
        count_matrix[i, i] = 0.0
    end
    normalization = sum(count_matrix, dims=1)
    holding_scale = 1 ./ mean.(holding_times(partitions, number_of_states; dt=dt))
    # calculate generator and handle edge case where no transitions occur
    for i in ProgressBar(eachindex(normalization))
        for j in generator_matrix[:, i].nzind
            generator_matrix[j, i] /= normalization[i]
        end
        generator_matrix[i, i] = -1.0
        for j in generator_matrix[:, i].nzind
            generator_matrix[j, i] *= holding_scale[i]
        end
        if normalization[i] == 0.0
            for j in generator_matrix[:, i].nzind
                generator_matrix[j, i] *= false
            end
        end
    end
    return generator_matrix
end
