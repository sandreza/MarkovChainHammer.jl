# two decompositions, ER and D, U 
using LinearAlgebra
export exit_rate, decomposition
"""
    exit_rate(Q)

# Description
    Decompose Generator Matrix Q into exit probabilities, E, and rates, R, such that Q = ER

# Arguments
- `Q::Matrix`: The generator matrix or transition matrix.

# Returns
- (; exit_probabilities, rates)::NamedTuple: The exit probabilities and rates of the generator matrix.
"""
function exit_rate(Q)
    rates = Diagonal([-Q[i,i] for i in 1:size(Q)[1]])
    exit_probabilities = Q * inv(Diagonal(rates))
    return (; exit_probabilities, rates)
end

"""
decomposition(Q)

# Description
    Decompose Generator Matrix Q into a reverible and volume preserving component. 

# Arguments
- `Q::Matrix`: The generator matrix or transition matrix.

# Returns
- (; exit_probabilities, rates)::NamedTuple: The exit probabilities and rates of the generator matrix.
"""
function decomposition(Q::AbstractArray)    
    _, V = eigen(Q)
    p = real.(V[:, end] ./ sum(V[:, end]))
    P = Diagonal(p)
    Qˢ = (Q + P * Q' * inv(P))/2
    Qᴬ = (Q - P * Q' * inv(P))/2
    return (; reversible = Qˢ, volume_preserving = Qᴬ)
end

