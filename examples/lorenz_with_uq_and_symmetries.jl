using ProgressBars, LinearAlgebra, Statistics
using MarkovChainHammer.BayesianMatrix
using MarkovChainHammer.TransitionMatrix: generator
using Random
Random.seed!(12345)
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

fixed_points = [[-sqrt(72), -sqrt(72), 27], [0.0, 0.0, 0.0], [sqrt(72), sqrt(72), 27]]
markov_states = fixed_points

timeseries = Vector{Float64}[]
markov_chain1 = Int64[]
initial_condition = [14.0, 20.0, 27.0]
push!(timeseries, initial_condition)
dt = 1.5586522107162 / 64
iterations = 64 * 16 * 10
for i in ProgressBar(2:iterations)
    local state = rk4(lorenz!, timeseries[i-1], dt)
    push!(timeseries, state)
    # partition state space according to most similar markov state
    markov_index = argmin([norm(state - markov_state) for markov_state in markov_states])
    push!(markov_chain1, markov_index)
end

# construct the Bayesian generator with an uninformative prior by default
Q1 = BayesianGenerator(markov_chain1; dt=dt)

# check against the empirical generator
Q = generator(markov_chain1; dt =dt)

println("The mean of our Bayesian generator is ")
display(mean(Q1))

println("The difference between the empirical and bayesian is ") 
display(mean(Q1) - Q)

println("however we now have uncertainties in the generator which we can check from the standard deviation")
display(std(Q1))

# Symmetries and Bayesian Update

function lorenz_symmetry(i)
    if i == 1
        return 3
    elseif i == 2
        return 2
    else
        return 1
    end
end

# Now we learn the symmetrized data
markov_chain2 = lorenz_symmetry.(markov_chain1)
Q2 = BayesianGenerator(markov_chain2; dt=dt)

# We take the posterior of the first simulation for the second simulation  
Q12 = BayesianGenerator(markov_chain2, Q1.posterior; dt = dt)
# and visa versa
Q21 = BayesianGenerator(markov_chain1, Q2.posterior; dt = dt)

# The two answers should be exactly equivalent
println("   ")
println("We can check that the order in which we assimilate data doesn't matter")
# and should correspond to learning all the data at once
println("The difference in means are  ") 
display(mean(Q12) - mean(Q21))
println("The difference in variances are  ") 
display(var(Q21) - var(Q12))


Q = Q12

# We can check the uncertainty on the eigenvalues
ΛVs = eigen_distribution(Q)
Λ1s = [ΛV.values[2] for ΛV in ΛVs] # the second eigenvalue

Λ̅ = mean(Λ1s)
Λ_of_mean =  eigvals(Q)[2]

println("We get different answers for the mean eigenvalue")
display(Λ̅)
println("and the eigenvalue of the mean")
display(Λ_of_mean)

## 
println("Although a bit of a heuristic we can check to see if a matrix constructed with more data is within the uncertainty of a matrix constructed with less data")

prior = MarkovChainHammer.BayesianMatrix.uninformative_prior(3)
Q_small_data = BayesianGenerator(markov_chain1[1:2000]; dt=dt)
Q_large_data = BayesianGenerator(markov_chain1; dt=dt)
Q_large_data = BayesianGenerator(markov_chain2, Q_large_data.posterior; dt=dt) 
# 

println("The mean of the small data generator is ")
display(mean(Q_small_data) )
println("The mean of the large data generator is ")
display(mean(Q_large_data) )
println("The 2 standard deviations of the small data matrix is ")
display(2 * std(Q_small_data) )
println("We see that the difference is with 2σ of the small data matrix")
display(abs.(mean(Q_small_data) - mean(Q_large_data)) .< 2*std(Q_small_data))
