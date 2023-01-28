import MarkovChainHammer.TransitionMatrix: holding_times, count_operator
using ProgressBars

struct ContinuousTimeEmpiricalProcess{S,T}
    holding_times::S
    cumulative_distribution::T 
end

"""
ContinuousTimeEmpiricalProcess(markov_chain; number_of_states=length(union(markov_chain)))

Construct a continuous time empirical process from a discrete time Markov chain.
The holding times are taken from the empirical distribution and the transition probabilities from the empirical transition probabilities

# Arguments 
- `markov_chain::Vector{Int}`: The discrete time Markov chain to construct the continuous time empirical process from.

# Keyword Arguments
- `number_of_states::Int=length(union(markov_chain))`: The number of states in the Markov chain.

# Returns
- `ContinuousTimeEmpiricalProcess`: A continuous time empirical process constructed from the discrete time Markov chain.

"""
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