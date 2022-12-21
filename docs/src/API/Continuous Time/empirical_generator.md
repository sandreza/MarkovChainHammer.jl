# [Empirical Generator](@id sec:empirical_generator)

In this section, we will see how to construct a generator from data. This is useful when we don't know the underlying generator. There is an inherent limitation in constructing generators that wasn't present when constructing the transfer operator: The data comes at a fixed timescale, thus the generator will only be known up to the timescales present in the timeseries. Nonetheless it is possible to construct a generator from the timeseries under a few assumptions. 

For example, suppose that we have the Markov chain with three states

```@example datadriven2
markov_chain = [1 1 2 2 2 1 3 1 1 2 2]
```

and we want to estimate the transfer operator from this data. We can do this by using the function `generator` which is imported from the `MarkovChainHammer.TransitionMatrix` module. This function takes a Markov chain as input and returns the empirical transition matrix of the Markov chain. We assume that the time intervals are uniform in time with timestep ``dt = 1``

```@example datadriven2
using MarkovChainHammer.TransitionMatrix: generator
dt = 1.0
Q = generator(markov_chain; dt = dt)
Q
```

This matrix is constructed by counting the amount of time spent in a state to construct the diagonal entries of the matrix, i.e., the holding times. We can manually import the holding times to understand the construction of the ``Q`` matrix, from the `MarkovChainHammer.TransitionMatrix` function `holding_times`

```@example datadriven2 
using MarkovChainHammer.TransitionMatrix: holding_times
state_holding_times = holding_times(markov_chain, 3; dt = dt)
```
The first argument to the function is the markov_chain, the second argument is the total number of states, and the keyword argument is the time step size associated with the entries of the markov_chain.

We see, from the markov chain, that state 1 spends 2 time units in state 1, followed by 1 time unit, then follow by 2 time units. State 2 spends 3 time units repeating itself, then two time units. And lastly, state 3 is only observed once for one time unit.

The diagonals of the ``Q`` matrix are given by taking the average of the holding times, then taking the reciprocal. We can verify this manually by importing the `mean` function from the `Statistics` package.

```@example datadriven2 
using Statistics
1 ./ mean.(state_holding_times)
```

We can verify the off diagonal terms as well. We just need to track the number of times that a given state was exited. For example, when leaving state 1 we observe ``1 \rightarrow 2``, read as 1 transitions to 2, ``1 \rightarrow 3``, and ``1 \rightarrow 2``. Thus going from state state 1 to state 2 is twice as likely as going from state 1 to state 3 and the associated probabilities are ``2/3`` for transitioning from 1 to 2 and ``1/3`` for transitioning from 2 to 3. These two probabilities then get multiplied by the reciprocal empircal holding time for being in state 1, which in this case is ``0.6`` to yield the first column of the ``Q`` matrix 

```@example datadriven2 
Q[:, 1]
```

For the other two states we observe that we only see state 2 transitioning to state 1 and state 3 only ever transitions to state 1, thus the matrix entries are ``Q_{32} = Q_{23} = 0``. The requirement that the columns must sum to zero then determines the other entries (or noticing that the probability of going to state 1 after having left state 2 is 1 and similarly for state 3).

This brings us to our first warning about constructing the empirical generator. It is most meaningful when there is enough timeseries data to stay within a given state for more than one "timestep". 

Our last comment is that changing ``dt`` amounts to rescaling the generator. That is to say, 

```@example datadriven2 
Q1 = generator(markov_chain; dt = 1.0)
Q2  = generator(markov_chain; dt = 2.0)
Q1 - 2.0 * Q2 
```

As a last comment we mention that one can also construct the transfer operator, then take the matrix logarithm, and then divide by ``dt`` to get another estimate of the generator; however, the resulting matrix no longer has an "infinitesimal" probabilistic interpretation, as the following example

```@example datadriven2 
using MarkovChainHammer.TransitionMatrix: perron_frobenius
dt = 1.0
ℳ = perron_frobenius(markov_chain)
log(ℳ) / dt
```

The columns still sum to zero, but the entries of the matrix no longer have an interpretation as holding times or probabilities. This example also shows why generators have a much more limited and special structure as compared to the transfer operator. The requirement that one generates probabilities over infinitesimal steps imposes severe restrictions. 