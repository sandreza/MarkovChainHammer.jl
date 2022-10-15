using ProgressBars, LinearAlgebra, Statistics
import MarkovChainHammer.TransitionMatrix: generator

# generate data
function lorenz!(ṡ, s)
    ṡ[1] = 10.0 * (s[2] - s[1])
    ṡ[2] = s[1] * (28.0 - s[3]) - s[2]
    ṡ[3] = s[1] * s[2] - (8 / 3) * s[3]
    return nothing
end

function rk4(f, s, dt)
    ls = length(s)
    k1 = zeros(ls)
    k2 = zeros(ls)
    k3 = zeros(ls)
    k4 = zeros(ls)
    f(k1, s)
    f(k2, s + k1 * dt / 2)
    f(k3, s + k2 * dt / 2)
    f(k4, s + k3 * dt)
    return s + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
end

timeseries = Vector{Float64}[]
markov_chain = Int64[]
s = [14.0, 20.0, 27.0]
subdiv = 200
push!(timeseries, s)
dt = 1.5586522107162 / subdiv
markov_states = [[-sqrt(72), -sqrt(72), 27], [0.0, 0.0, 0.0], [sqrt(72), sqrt(72), 27]]
iterations = 1000000
for i in ProgressBar(2:iterations)
    local s = rk4(lorenz!, timeseries[i-1], dt)
    push!(markov_chain, argmin([norm(s - ms) for ms in markov_states]))
    push!(timeseries, s)
end

## construct transition matrix
pQ = generator(markov_chain; dt=dt)
Q = copy(pQ)
# symmetrize
Q[1, 1] = Q[3, 3] = (pQ[1, 1] + pQ[3, 3]) * 0.5
Q[2, 1] = Q[2, 3] = (pQ[2, 1] + pQ[2, 3]) * 0.5
Q[3, 1] = Q[1, 3] = (pQ[3, 1] + pQ[1, 3]) * 0.5
Q[1, 2] = Q[3, 2] = (pQ[1, 2] + pQ[3, 2]) * 0.5
# pick out the steady state distribution
Λ, V = eigen(Q)
p = V[:, end] ./ sum(V[:, end])

## look ensemble average versus time average of an observable
observable(u) = u[3]
# ensemble average
g_ensemble =  sum(observable.(markov_states) .* p)
g_timeseries = mean(observable.(timeseries))
println("The ensemble average is $(g_ensemble)")
println("The timeseries average is $(g_timeseries)")

