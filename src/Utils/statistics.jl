using LinearAlgebra, Statistics, ProgressBars

export steady_state, scaled_entropy, koopman_modes, decorrelation_times, autocovariance


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
`scaled_entropy(p)`

# Description
    Calculate the scaled_entropy of a probability distribution. Normalizethe entropy of the empirical distribution to that of the uniform distribution.

# Arguments
- `p::Vector`: The probability distribution.

# Returns
- `entropy_value::Real`: The entropy of the probability distribution.
"""
function scaled_entropy(p)
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
    Λ, V = eigen(Q')
    return V
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

"""
    autocovariance(x; timesteps=length(x), progress = false)

Calculate the autocovariance of a timeseries `x` with `timesteps` lags.

# Arguments
- `x::AbstractVector`: timeseries

# Keyword Arguments
- `timesteps::Int`: number of lags
- `progress::Bool`: show a progress bar

# Returns
- `autocov::Vector`: autocovariance of `x` with `timesteps` lags
"""
function autocovariance(x; timesteps=length(x), progress = false)
    μ = mean(x)
    autocor = zeros(timesteps)
    progress ? iter = ProgressBar(1:timesteps) : iter = 1:timesteps
    for i in iter
        autocor[i] = mean(x[i:end] .* x[1:end-i+1]) - μ^2
    end
    return autocor
end

"""
    autocovariance(g⃗, Q::Eigen, timelist; progress=false)

Calculate the autocovariance of observable g⃗ with generator matrix Q and times at timelist. 

# Arguments
- `g⃗::AbstractVector`: observable
- `Q::Eigen`: eigenvalue decomposition of generator matrix
- `timelist::AbstractVector`: times at which to calculate autocovariance

# Keyword Arguments
- `progress::Bool=false`: show a progress bar

# Returns
- `autocov::Vector`: autocovariance of observable g⃗ with generator matrix Q and times at timelist
"""
function autocovariance(g⃗, Q::Eigen, timelist; progress=false)
    @assert all(real.(Q.values[1:end-1]) .< 0) "Did not pass an ergodic generator matrix"

    autocov = zeros(length(timelist))
    # Q  = V Λ V⁻¹
    Λ, V = Q
    p = real.(V[:, end] ./ sum(V[:, end]))
    v1 = V \ (p .* g⃗)
    w1 = g⃗' * V
    μ = sum(p .* g⃗)
    progress ? iter = ProgressBar(eachindex(timelist)) : iter = eachindex(timelist)
    for i in iter
        autocov[i] = real(w1 * (exp.(Λ .* timelist[i]) .* v1) - μ^2)
    end
    return autocov
end


"""
    autocovariance(g⃗, Ps::Vector{Matrix{Float64}}, steps::Int; progress = false)

Calculate the autocovariance of observable g⃗ with transition matrices Ps and number of steps.

# Arguments
- `g⃗::AbstractVector`: observable
- `Ps::Vector{Matrix{Float64}}`: transition matrices
- `steps::Int`: number of steps

# Keyword Arguments
- `progress::Bool=false`: show a progress bar

# Returns
- `autocov::Vector`: autocovariance of observable g⃗ with transition matrices Ps and number of steps
"""
function autocovariance(observable, Ps::Vector{Matrix{Float64}}, steps::Int; progress = false)
    @assert length(Ps) >= steps "The length of Ps should be greater than or equal to steps"
    autocor = zeros(steps + 1)
    p = steady_state(Ps[1])
    μ² = sum(observable .* p)^2
    autocor[1] = observable' * (observable .* p) - μ²
    progress ? iter = ProgressBar(1:steps) : iter = 1:steps
    for i in iter
        autocor[i+1] = observable' * Ps[i] * (observable .* p) - μ²
    end
    return autocor
end

autocovariance(g⃗, Q, timelist; progress = false) = autocovariance(g⃗, eigen(Q), timelist; progress = progress)
