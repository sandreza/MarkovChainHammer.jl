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

function generate(matrix, n, dt, initial_condition)
    if all(sum(matrix, dims = 1) .≈ 1)
        P = matrix 
    else
        P = exp(matrix * dt)
    end
    cP = cumsum(P, dims=1)
    markov_chain = zeros(Int, n)
    markov_chain[1] = initial_condition
    for i = 2:n
        markov_chain[i] = next_state(markov_chain[i-1], cP)
    end
    return markov_chain
end

generate(Q, n; dt=1, initial_condition=rand(1:size(Q)[1])) = generate(Q, n, dt, initial_condition)

function generate(Q; dt = 1, decorrelation_threshold = 0.01, initial_condition=rand(1:size(Q)[1]) )
    n = step_heuristic(Q; dt = 1, decorrelation_threshold = decorrelation_threshold)
    return generate(Q, n, dt, initial_condition)
end

function step_heuristic(matrix; dt = 1, decorrelation_threshold = 0.01)
    if all(sum(matrix, dims = 1) .≈ 1)
        P = matrix 
    else
        P = exp(matrix * dt)
    end
    ll = eigvals(P)
    real_ll = abs.(real.(ll))
    
    # edgecase
    if length(real_ll[real_ll .< 1 - 100* eps()]) == 0
        # since not ergodic anyways
        return 2
    end
    slowest = maximum(real_ll[real_ll .< 1 - 100* eps()])
    n = ceil(Int, log(decorrelation_threshold) / log(slowest)) * 10000
    if n < 10^4 
        return 10^4 
    elseif n > 10^7 
        println("warning: step_heuristic is returning a large number of steps.")
        return 10^7 
    else
        return n
    end
end
