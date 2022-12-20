# [Transfer Operators and Markov Chains](@id sec:transfer_operators)

Markov chains are a stochastic process whose future state only depends on the current state. In this repository we only consider Markov chains with a finite state space, thus the transition probabilities are characterized by matrix. In this section we will see how to generate a Markov chain from a known transition matrix. The transition matrix is also known as the transfer operator or the Perron-Frobenius matrix / operator. It can also be viewed as the adjoint of the [Koopman operator](https://en.wikipedia.org/wiki/Composition_operator).

The **convention** taken in this repository is that all transfer operators are column stochastic. For example, the following ``2 \times 2``  [column stochastic matrix](https://en.wikipedia.org/wiki/Stochastic_matrix) characterizes a Markov chain made up of 2 discrete states,

```math
\begin{aligned}
    \mathcal{M} =
    \begin{bmatrix}
    0.6 & 0.3 \\
    0.4 & 0.7
    \end{bmatrix}.
\end{aligned}
```

The first column, ``c_1``,

```math
\begin{aligned}
    c_1 =
    \begin{bmatrix}
    0.6 \\
    0.4 
    \end{bmatrix},
\end{aligned}
```

gives information about what happens to a link in the chain that is state 1. The probability, ``\mathcal{M}_{11}``, of staying in state 1 given that we are in state 1 is ``\mathcal{M}_{11} = 0.6`` and the probability of going to state 2 given that we are in state 1 is ``\mathcal{M}_{21} = 0.4``. Similarly, the second column, ``c_2``,

```math
\begin{aligned}
    c_2 =
    \begin{bmatrix}
    0.3 \\
    0.7
    \end{bmatrix},
\end{aligned}
```

gives information about what happens to a link in the chain that is state 2. The probability, ``\mathcal{M}_{22}``, of staying in state 2 given that we are in state 2 is ``\mathcal{M}_{22} = 0.7`` and the probability of going to state 1 given that we are in state 2 is ``\mathcal{M}_{12} = 0.3``. 

A sample markov chain is constructed from a transfer operator using *generate* function from the *MarkovChainHammer.Trajectory* module. The following code snippet constructs 10 steps of a Markov chain from the transfer operator ``\mathcal{M}`` above.

```@example generate_chain
using MarkovChainHammer.Trajectory: generate
using Random #hide
Random.seed!(1234) #hide
ℳ = [0.6 0.3; 0.4 0.7]
steps = 10
markov_chain = generate(ℳ, steps)'
```

The Markov chain is but one possible realization of the stochastic process. If we were to run the generate function again, we would get a different realization of the Markov chain,

```@example generate_chain
markov_chain = generate(ℳ, steps)'
```

If the number of steps is not specified, then the code attempts to generated a chain with roughly 10,000 independent samples of the process based on a decorrelation time threshold,  


```@example generate_chain
markov_chain = generate(ℳ)'
```