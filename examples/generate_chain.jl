using LinearAlgebra: eigen, norm
using Random

using MarkovChainHammer
using MarkovChainHammer.Trajectory: generate
using MarkovChainHammer.TransitionMatrix: generator, perron_frobenius
using MarkovChainHammer.TransitionMatrix: steady_state 

Random.seed!(123456789)

# Markov Chain Transition Matrices
Q = [-1/2 1/3 0.0; 1/2 -2/3 1/3; 0.0 1/3 -1/3]
dt = 0.1
P = exp(Q * dt)
p_exact = steady_state(P)

# Generate Markov Chain
markov_chain = generate(P)

# Generate Empirical Construction from Markov Chain
p_empirical = [sum(markov_chain .== i) / length(markov_chain) for i in 1:3]
Q_empirical = generator(markov_chain; dt=0.1)
P_empirical = perron_frobenius(markov_chain)

# Compare
println("Exact steady state: "); 
display(p_exact)
println("Empirical steady state: "); 
display(p_empirical)
println("Relative Percent Error: ", 100 * norm(p_exact - p_empirical) / norm(p_exact))
println("--------------------------------------")

# Compare
println("Exact transfer operator: "); 
display(P)
println("Empirical transfer operator: "); 
display(P_empirical)
println("Relative Percent Error: ", 100 * norm(P - P_empirical) / norm(P))
println("--------------------------------------")

##
# Compare
println("Exact generator: "); 
display(Q)
println("Empirical generator: "); 
display(Q_empirical)
println("Relative Percent Error: ", 100 * norm(Q - Q_empirical) / norm(Q))
println("--------------------------------------")