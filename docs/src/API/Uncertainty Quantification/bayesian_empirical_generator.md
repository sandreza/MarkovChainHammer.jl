# [Bayesian Empirical Generator](@id sec:bayesian_empirical)

In this section we will show how to utilise the Bayesian empirical generator to estimate the uncertainty of the entries in the empirical generator from finite data. We we will do so by generating a Markov chain from a generator matrix. First let's load a few packages and functions that we've become familiar with over the last few sections.
```@example bayesian_empirical_generator
using Random, MarkovChainHammer, Statistics
using MarkovChainHammer.TransitionMatrix: generator
using MarkovChainHammer.Trajectory: generate
Random.seed!(1234)
```
Observe the use of the Random package to set the seed for the random number generator. This is important for reproducibility of the results. A seed should only be defined once at the beginning of a script. 

We will now create a generator matrix Q and generate a markov process.

```@example bayesian_empirical_generator
Q = [-1.0 4/3 2; 1/4 -2.0 1; 3/4 2/3 -3.0]
dt = 0.01
markov_chain = generate(Q, 10000; dt = dt)'
```

We first construct the empirical generator as before to serve as a comparison point for the Bayesian empirical generator which will come later.

```@example bayesian_empirical_generator
Qempirical = generator(markov_chain; dt = dt)
```

We are now ready to introduce the Bayesian Generator matrix. We first load in the BayesianMatrix submodule.

```@example bayesian_empirical_generator
using MarkovChainHammer.BayesianMatrix
```

This submodule naturally exports a number of structs and functions. We can see the names of these functions by typing the following: 
```@example bayesian_empirical_generator
names(MarkovChainHammer.BayesianMatrix)
```

In particular a BayesianGenerator object is exported. We use the BayesianGenerator just like a normal generator in order to construct a BayesianGenerator

```@example bayesian_empirical_generator
Q_bayes = BayesianGenerator(markov_chain; dt = dt)
```

This is no longer a regular matrix, but rather a random matrix whose entries are given by a probability distribution consistent with finite sampling from a markov process. In the present context (this will change later when we discuss prior distributions) the mean of the Bayesian matrix reproduces the same matrix as the emperically obtained matrix, as we can check

```@example bayesian_empirical_generator
mean(Q_bayes) - Qempirical
```

However, we have more than just a mean value for the Bayesian matrix, we also have variances 
```@example bayesian_empirical_generator
var(Q_bayes)
```
and standard deviations
```@example bayesian_empirical_generator
std(Q_bayes)
```

We can check that the mean falls within two standard deviations of the true empirical estimate with high probability

```@example bayesian_empirical_generator
abs.(mean(Q_bayes) - Q) .<  2 * std(Q_bayes)
```

