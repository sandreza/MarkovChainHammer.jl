# [Bayesian Empirical Generator](@id sec:bayesian_empirical)

```@example bayesian_empirical_generator
using Random, MarkovChainHammer
using MarkovChainHammer.TransitionMatrix: generator
using MarkovChainHammer.Trajectory: generate
Random.seed!(1234)
Q = [-1.0 4/3 2; 1/4 -2.0 1; 3/4 2/3 -3.0]
dt = 0.01
markov_chain = generate(Q, 10000; dt = dt)
Qempirical = generator(markov_chain; dt = dt)
```


```@example bayesian_empirical_generator
using MarkovChainHammer.BayesianMatrix
```

```@example bayesian_empirical_generator
names(MarkovChainHammer.BayesianMatrix)
```

```@example bayesian_empirical_generator
Q_bayes = BayesianGenerator(markov_chain; dt = dt)
```
