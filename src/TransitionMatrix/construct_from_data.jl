using Statistics

function count_operator(markov_chain, number_of_states)
    count_matrix = zeros(number_of_states, number_of_states)
    for i in 1:length(markov_chain)-1
        count_matrix[markov_chain[i+1], markov_chain[i]] += 1
    end
    return count_matrix
end

count_operator(markov_chain) = count_operator(markov_chain, maximum(markov_chain))

function perron_frobenius(markov_chain, number_of_states)
    count_matrix = count_operator(markov_chain, number_of_states)
    normalization = sum(count_matrix, dims=1)
    perron_frobenius_matrix = count_matrix ./ normalization
    # handle edge case where no transitions occur
    for i in eachindex(normalization)
        if normalization[i] == 0.0
            perron_frobenius_matrix[:, i] *= false
            perron_frobenius_matrix[i, i] = 1.0
        end
    end
    return perron_frobenius_matrix
end

perron_frobenius(markov_chain) = perron_frobenius(markov_chain, maximum(markov_chain))

function holding_times(markov_chain, number_of_states; dt=1)
    holding_times = [[] for n in 1:number_of_states]
    push!(holding_times[markov_chain[1]], dt)
    M = length(markov_chain)
    for i in 2:M
        current_state = markov_chain[i]
        previous_state = markov_chain[i-1]
        if current_state == previous_state
            holding_times[current_state][end] += dt
        else
            push!(holding_times[current_state], dt)
        end
    end
    return holding_times
end

holding_times(markov_chain; dt=1) = holding_times(markov_chain, maximum(markov_chain); dt=dt)

function generator(markov_chain, number_of_states; dt=1)
    # calculate average time spent in a state
    holding_time = holding_times(markov_chain, number_of_states, dt=dt)
    holding_scale = zeros(length(holding_time))
    for i in eachindex(holding_scale)
        if length(holding_time[i]) > 0
            holding_scale[i] = 1 / mean(holding_time[i])
        end
    end
    # calculate transitions away from state
    count_matrix = count_operator(markov_chain, number_of_states)
    for i in 1:number_of_states
        count_matrix[i, i] = 0.0
    end
    normalization = sum(count_matrix, dims=1)
    # calculate generator and handle edge case where no transitions occur
    generator_matrix = count_matrix ./ normalization
    for i in eachindex(normalization)
        generator_matrix[i, i] = -1.0
        generator_matrix[:, i] *= holding_scale[i]
        if normalization[i] == 0.0
            generator_matrix[:, i] *= false
        end
    end
    return generator_matrix
end

generator(markov_chain; dt=1) = generator(markov_chain, maximum(markov_chain); dt=dt)

function symmetric_generator(markov_chain, symmetries, number_of_states; dt=1)
    Q = generator(markov_chain, number_of_states; dt=dt)
    # grab holding times
    diagQ⁻¹ = 1 ./ reshape([Q[i, i] for i in 1:number_of_states], (1, number_of_states))
    # normalize so that we have probabilities
    Q .*= diagQ⁻¹
    for symmetry in symmetries
        tmpQ = generator(symmetry.(markov_chain), number_of_states; dt=dt)
        tmp_diagQ⁻¹ = 1 ./ reshape([tmpQ[i, i] for i in 1:number_of_states], (1, number_of_states))
        Q .+= (tmpQ .* tmp_diagQ⁻¹)
        # need to recompute average holding time
        diagQ⁻¹ .+= tmp_diagQ⁻¹
    end
    # average over symmetries 
    # note, there is a factor 1 + length(symmetries) that cancels out
    Q ./= diagQ⁻¹
    return Q
end

symmetric_generator(markov_chain, symmetries; dt=1) = symmetric_generator(markov_chain, symmetries, length(union(markov_chain)); dt=dt)

function symmetric_perron_frobenius(markov_chain, symmetries, number_of_states)
    P = perron_frobenius(markov_chain, number_of_states)
    for symmetry in symmetries
        P .+= perron_frobenius(symmetry.(markov_chain), number_of_states)
    end
    P .*= 1 / (1 + length(symmetries))
    return P
end

symmetric_perron_frobenius(markov_chain, symmetries) = symmetric_perron_frobenius(markov_chain, symmetries, length(union(markov_chain)))
