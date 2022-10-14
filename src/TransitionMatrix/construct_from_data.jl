using Statistics

function holding_times(markov_chain, number_of_states; dt=1)
    holding_times = [[] for n in 1:number_of_states]

    push!(holding_times[markov_chain[1]], dt)

    M = length(markov_chain)
    # loop through steps
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
    # dynamic container for holding times
    holding_times = [[] for n in 1:number_of_states]
    # dynamic container for states
    state_list = []

    push!(holding_times[markov_chain[1]], dt)
    push!(state_list, markov_chain[1])

    M = length(markov_chain)
    # loop through steps
    for i in 2:M
        current_state = markov_chain[i]
        previous_state = markov_chain[i-1]
        if current_state == previous_state
            holding_times[current_state][end] += dt
        else
            push!(state_list, current_state)
            push!(holding_times[current_state], dt)
        end
    end

    # construct transition matrix from state list 
    constructed_T = zeros(number_of_states, number_of_states)

    # first count how often state=i transitions to state=j 
    # the current state corresponds to a column of the matrix 
    # the next state corresponds to the row of the matrix
    number_of_unique_states = length(state_list)
    for i in 1:number_of_unique_states-1
        local current_state = state_list[i]
        local next_state = state_list[i+1]
        constructed_T[next_state, current_state] += 1
    end

    # edge case. Only see one state
    if all(constructed_T .== 0.0)
        return constructed_T
    end

    # now we need to normalize
    normalization = sum(constructed_T, dims=1)
    normalized_T = constructed_T ./ normalization

    # now we account for holding times 
    holding_scale = zeros(length(holding_times))
    for i in eachindex(holding_scale)
        if length(holding_times[i]) > 0
            holding_scale[i] = 1 / mean(holding_times[i])
        end
    end

    # need to handle edge case where no transitions occur
    for i in eachindex(holding_scale)
        if normalization[i] == 0.0
            normalized_T[:, i] *= false
        else
            normalized_T[i, i] = -1.0
            normalized_T[:, i] *= holding_scale[i]
        end
    end
    return normalized_T
end

generator(markov_chain; dt=1) = generator(markov_chain, maximum(markov_chain); dt=dt)


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
        end
    end
    return perron_frobenius_matrix
end

perron_frobenius(markov_chain) = perron_frobenius(markov_chain, maximum(markov_chain))
