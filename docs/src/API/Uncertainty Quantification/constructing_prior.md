# [Constructing Prior Distributions](@id sec:constructing_priors)

In the previous section we saw how to use the BayesianGenerator object, similar to how we use the generator function. We now introduce the GeneratorParameterDistribution object in order to construct prior distributions for Bayesian inference. 

We first load the necessary packages and functions

```@example constructing_prior
using Random, MarkovChainHammer, Statistics
using MarkovChainHammer.BayesianMatrix
using MarkovChainHammer.TransitionMatrix: generator
using MarkovChainHammer.Trajectory: generate
Random.seed!(1234)
```

We now introduce the GeneratorParameterDistribution object, which was exported from MarkovChainHammer.BayesianMatrix

```@example constructing_prior
number_of_states = 3
prior = GeneratorParameterDistributions(number_of_states)
```

The prior distribution constitutes a guess for the entries of the matrix. We can check the mean of a random matrix generator from the prior distribution by calling the mean function on the prior
```@example constructing_prior
mean(prior)
```

From a generator we now create a markov chain

```@example constructing_prior
Q = [-1.0 4/3 2; 1/4 -2.0 1; 3/4 2/3 -3.0]
dt = 0.01
markov_chain = generate(Q, 10000; dt = dt)'
```

and now use our prior distribution along with the BayesianGenerator object

```@example constructing_prior
Q_bayes = BayesianGenerator(markov_chain, prior; dt = dt)
```

We see that we get a different answer than what we had before due to the presence of the prior distribution

```@example constructing_prior
Q_bayes_uninformative_prior = BayesianGenerator(markov_chain; dt = dt)
```

In particular we can check the mean 
```@example constructing_prior
mean(Q_bayes) - mean(Q_bayes_uninformative_prior)
```
