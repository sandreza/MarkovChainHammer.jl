# NEED ERROR HANDLING
import MarkovChainHammer.Trajectory: next_state, generate

struct ContinuousTimeEmpiricalProcess{S,T}
    holding_times::S
    cumulative_distribution::T # change to cumulative sum for simulation purposes
end

function ContinuousTimeEmpiricalProcess(markov_chain; number_of_states=length(union(markov_chain)))
    ht = holding_times(markov_chain)
    count_matrix = count_operator(markov_chain, number_of_states)
    count_matrix = count_matrix - Diagonal(count_matrix)
    Ntotes = sum(count_matrix, dims=1)
    if any(Ntotes .== 0)
        throw(ArgumentError("Some states are not reachable from the initial condition."))
    end
    pmatrix = count_matrix ./ Ntotes
    cP = cumsum(pmatrix, dims=1)
    return ContinuousTimeEmpiricalProcess(ht, cP)
end

function generate(process::ContinuousTimeEmpiricalProcess, n, initial_condition)
    simulated_chain = Int64[]
    push!(simulated_chain, initial_condition)
    for i in ProgressBar(1:n)
        current_state = simulated_chain[end]
        htempirical = rand(process.holding_times[current_state])
        for i in 1:htempirical
            push!(simulated_chain, current_state)
        end
        simulated_chain = vcat(simulated_chain...)
        push!(simulated_chain, next_state(current_state, process.cumulative_distribution))
    end
    return simulated_chain
end