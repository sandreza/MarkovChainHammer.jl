# [Data-driven Transfer Operators](@id sec:empirical_transfer_operators)

[Previously](@ref sec:transfer_operators), we have seen how to construct a transfer operator from a Markov chain. In this section, we will see how to construct a transfer operator from data. This is useful when we don't know the underlying transfer operator.

For example, suppose that we have the Markov chain

```@example datadriven
markov_chain = [1 1 2 2 2 1]
```

and we want to estimate the transfer operator from this data. We can do this by using the function `perron_frobenius` which is imported from the `MarkovChainHammer.TransitionMatrix` module. This function takes a Markov chain as input and returns the empirical transition matrix of the Markov chain.

```@example datadriven
using MarkovChainHammer.TransitionMatrix: perron_frobenius
P = perron_frobenius(markov_chain)
P
```

This matrix is constructed by counting the number of times that the chain transitions from one state to another and dividing by the total number of transitions. For example, the entry `P[1, 1]` is the number of times that the chain transitions from state 1 to state 1 divided by the total number of transitions away from state 1. In this case, the chain transitions from state 1 to state 1, a total of 1 time. Furthermore, the chain transitions from state 1  to state 2 a total of 1 time. Thus the total number of transitions is 2 and the entries of the first column are ```P[1, 1] = 1/2``` and ```P[2, 1] = 1/2```.

The second column of the matrix is constructed in the same way. The chain transitions from state 2 to state 1 a total of 1 time and from state 2 to state 2 a total of 2 times. Thus the entries of the second column are ```P[1, 2] = 1/3``` and ```P[2, 2] = 2/3```.
