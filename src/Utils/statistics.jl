using LinearAlgebra, Statistics, ProgressBars
import MarkovChainHammer.TransitionMatrix: steady_state

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
        autocov[i] = real(w1 * (exp.(Λ .* timelist[i]) .* v1)) - μ^2
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
