# [Overview](@id sec:api_overview)

This repository gathers various convenience tools for analyzing and generating [Markov Chains](https://en.wikipedia.org/wiki/Markov_chain) with a finite state space. The following subsections increase in conceptual complexity and are intended to be read in order.

## Discrete time Markov chains
The following sections contain a review of Markov Chains and an introduction to functionality in the repository. 
0. [Transfer Operators and Markov Chains](@ref sec:transfer_operators): Generate a Markov Chain from a transfer operator
0. [Data-driven Transfer Operators](@ref sec:empirical_transfer_operators): Construct empirical transition matrices and distributions from a Markov Chain
0. [Convergence of Transfer Operators](@ref sec:to_convergence): Verify convergence of empirical transition matrices and distributions to the true transition operator

## Continuous time Markov chains
The following sections contain more advanced mathematical concepts using continuous time Markov processes 
0. [Generators and Markov Chains](@ref sec:generators): Create a Markov chain from a generator matrix
0. [Data-driven Generator](@ref sec:empirical_generator): Construct empirical generators from a Markov Chain
0. [Convergence of Generators](@ref sec:generator_convergence): Verify convergence of empirical generators to the true generator


