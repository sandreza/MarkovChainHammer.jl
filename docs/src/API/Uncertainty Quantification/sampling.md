# [Sampling from the Bayesian Generator](@id sec:generator_sampling)

We now discuss how to use the BayesianGenerator to propogate uncertainty. 

As usual we start by bringing in packages

```@example sampling
using Random, MarkovChainHammer, Statistics, LinearAlgebra
using MarkovChainHammer.BayesianMatrix
using MarkovChainHammer.TransitionMatrix: generator
using MarkovChainHammer.Trajectory: generate
Random.seed!(1234)
```

and starting off with a generator that generates a stochastic process
```@example sampling
Q = [-1.0 4/3 2; 1/4 -2.0 1; 3/4 2/3 -3.0]
dt = 0.01
markov_chain = generate(Q, 10000; dt = dt)'
```

We construct the BayesianGenerator object from a prior distribution
```@example sampling
number_of_states = 3
prior = GeneratorParameterDistributions(number_of_states)
Q_bayes = BayesianGenerator(markov_chain, prior; dt = dt)
```

We can now sample from the BayesianGenerator object which gives a random matrix drawn from the posterior distribution
```@example sampling
rand(Q_bayes)
```
If we call the rand function again we get a different realization
```@example sampling
rand(Q_bayes)
```
We can call the rand function with an additional integer argument to get a list of realizations
```@example sampling
number_of_samples = 100
Q_bayes_list = rand(Q_bayes, number_of_samples)
```
from whence we can compute the sample mean and compare to the analytic mean
```@example sampling
mean(Q_bayes_list) - mean(Q_bayes)
```
and similarly for the sample variance
```@example sampling
var(Q_bayes_list) - var(Q_bayes)
```

