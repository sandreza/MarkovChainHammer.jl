using LinearAlgebra

"""
`steady_state(Q)`

# Description
    Calculate the steady state of a generator matrix.

# Arguments
- `Q::Matrix`: The generator matrix or transition matrix.

# Returns
- `p::Vector`: The steady state of the generator matrix.
"""
function steady_state(Q)
    Λ, V = eigen(Q)
    return real(V[:, end] ./ sum(V[:,end]))
end

"""
`entropy(p)`

# Description
    Calculate the entropy of a probability distribution. Normalized by the entropy of the uniform distribution.

# Arguments
- `p::Vector`: The probability distribution.

# Returns
- `entropy_value::Real`: The entropy of the probability distribution.
"""
function entropy(p)
    N = length(p)
    entropy_value = 0.0
    for i in 1:N
        if p[i] > 100 * eps(eltype(p))
            entropy_value += p[i] .* log.(p[i])
        end
    end
    return -entropy_value / log(N)
end

"""
`koopman_modes(Q)`

# Description
    Calculate the koopman modes of a generator matrix.

# Arguments
- `Q::Matrix`: The generator matrix or transition matrix.

# Returns
- `V⁻¹::Matrix`: The koopman modes of the generator matrix. Each row is a koopman mode
"""
function koopman_modes(Q)
    Λ, V = eigen(Q)
    V⁻¹ = inv(V)
    return V⁻¹
end

"""
`decorrelation_times(Q)`

# Description
    Calculate the decorrelation times of a generator matrix.

# Arguments
- `Q::Matrix`: The generator matrix or transition matrix.

# Returns
- `decorrelation_times::Vector`: The decorrelation times of the generator matrix.
"""
function decorrelation_times(Q)
    Λ = eigvals(Q)
    return 1 ./ real.(Λ)
end