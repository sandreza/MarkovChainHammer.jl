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

```@example sampling
Q = [-1.0 4/3 2; 1/4 -2.0 1; 3/4 2/3 -3.0]
dt = 0.01
markov_chain = generate(Q, 10000; dt = dt)'
```


```@example sampling
number_of_states = 3
prior = GeneratorParameterDistributions(number_of_states)
Q_bayes = BayesianGenerator(markov_chain, prior; dt = dt)
```


```@example sampling
rand(Q_bayes)
```

```@example sampling
rand(Q_bayes)
```

```@example sampling
number_of_samples = 100
Q_bayes_list = rand(Q_bayes, number_of_samples)
```

```@example sampling
mean(Q_bayes_list) - mean(Q_bayes)
```

```@example sampling
var(Q_bayes_list) - var(Q_bayes)
```

