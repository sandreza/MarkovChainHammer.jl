using Distributions, LinearAlgebra, Random

# Define the jump map
function next_state(current_state_index::Int, cT)
    vcT = view(cT, :, current_state_index)
    u = rand(Uniform(0, 1))
    # choose a random uniform variable and decide next state
    # depending on where one lies on the line with respect to the probability 
    # of being found between given probabilities
    for i in eachindex(vcT)
        if u < vcT[i]
            return i
        end
    end

    return i
end

function generate(Q, n; dt=1, initial_condition=rand(1:size(Q)[1]))
    expQ = exp(dt * Q) # get the transition probabilities
    cexpQ = cumsum(expQ, dims=1)
    markov_chain = zeros(Int, n)
    markov_chain[1] = initial_condition
    for i = 2:n
        markov_chain[i] = next_state(markov_chain[i-1], cexpQ)
    end
    return markov_chain
end
