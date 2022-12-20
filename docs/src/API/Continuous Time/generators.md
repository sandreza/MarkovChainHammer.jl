# [Generators](@id sec:generators) 

In this section we discuss continuous time markov processes and methods for generating Markov chains from them. They are similar to a discrete time Markov chain, but with a notion of time built into their construction. The idea is to break up the time interval into a sequence of discrete time steps and then construct a discrete time Markov chain from the continuous time Markov process.

We start by building from the knowledge of the discrete time Markov chain and [transfer operators](@ref sec:transfer_operators). As mentioned before the columns of the transfer operator give information on where to go at a subsequent discrete step. But now we want to shift and develop a method that is consistent with the notion of continuous time. Thus the first thing we need to think about is, "How far do we go in time?" and then secondly "How do we make this consistent with the notion of a Markov process?". 


## Motivation
To motivate the use of a generator, we reexamine the transfer operator

```@example generator
ℳ = [0.6 0.3; 0.4 0.7]
ℳ
```

The transfer operator determines the probability of state j (here j = 1 or 2) to move to state i in one step. But suppose that we want to know the probability of moving from one state to another in *two* steps. Suppose that we are in state ``s_0`` at time 0. For example, if we were in state 2, then 
```@example generator
s₀ = [1, 0]
```
Thus the probability of being in state 1 or state 2 at time 1 is given by matrix multiplication with the state 
```@example generator
s₁ = ℳ * s₀
```
Now we would find that 60% of our ensembles would be in state 1 after one step and 40% would be in state 2 after one step. But now suppose that we want to figure out how many of our ensembles of states would be in state 1 or state 2 after another step. We simply apply the operator again 
```@example generator
s₂ = ℳ * s₁
```
To get that 48% of our states ended up in state 1 and 52% of our states ended up in state 2. In total, if we want to just jump straight from our initial state ``s_0`` to state ``s_2``, we can just combine the two steps to get 
```@example generator
ℳ * ℳ * s₀
```

Thus we see if we want to move one step, we just multiply ``s_0`` by ``ℳ``, if we want to move two steps we multiply ``s_0`` by, ``\mathcal{ℳ}^2``. In general, if we want to move ``n``-steps, where ``n`` is a natural number, then we mulitply by ``\mathcal{ℳ}^n``. If we want to move "0" steps, we can use the convention that ``\mathcal{ℳ}^0 = \mathbb{I}`` where ``\mathbb{I}`` is the identity matrix.

But now suppose that we don't want to compute at just integer steps, we want to compute what "happens in between", that is to say, raise to a non-integer power. One way to define this for matrices is to make use of the identity
```math
\mathcal{M}^n = \exp[n \log(\mathcal{M})]
```
where the log and exp functions are the [matrix logarithm](https://en.wikipedia.org/wiki/Logarithm_of_a_matrix#:~:text=In%20mathematics%2C%20a%20logarithm%20of,function%20of%20the%20matrix%20exponential.) and [matrix exponential](https://en.wikipedia.org/wiki/Matrix_exponential). For example, the matrix logarithm of ``\mathcal{M}`` is
```@example generator
logℳ = log(ℳ)
logℳ = real.(log(ℳ)) #hide
``` 

and we can check 
```@example generator
exp(logℳ) - ℳ
real.(exp(real.(log(ℳ))) - ℳ) # hide
``` 
We see that the mathematical relation holds to [machine precision](https://en.wikipedia.org/wiki/Machine_epsilon).

We can now check the formula for ``n=2``
```@example generator
n = 2
exp(n * log(ℳ)) - ℳ^n
real.(exp(n * log(ℳ)) - ℳ^n) # hide
``` 
which we again see holds to [machine precision](https://en.wikipedia.org/wiki/Machine_epsilon).

And finally, we can see what happens if we take non-integer steps
```@example generator
n = 1.5
ℳhalf = exp(n * log(ℳ)) 
``` 

and check to see taking two steps of size 1.5 is the same as taking one step of size 3. We square the matrix corresponding to a step size of 1.5 and check that it is the same as the matrix corresponding to a step size of 3. We take the difference of the two calculations
```@example generator
ℳhalf^2 - ℳ^3
``` 
to see that the entries are the same to machine precision.

This brings us to the point we were can finally introduce the generator. Presently, the generator is the matrix logarithm of the transfer operator 

```@example generator
Q = log(ℳ)
Q = real.(log(ℳ)) # hide
``` 

However, we have actually presented the relation backwards! The transfer operator for a time ``\tau`` in the future is the matrix exponential of the generator Q

```@example generator
P(τ) = exp(Q * τ)
nothing # hide
``` 

In order for the formula to hold for all times ``\tau``, stricter requirements must be imposed on the generator. In order for ``P(\tau)`` to be a probability for all times ``\tau``, we require 
0. The sum of the columns of Q must be zero 
0. The off diagonal terms of Q must be greater than or equal to zero

These two requirements can be derived by taking ``\tau`` to be an infinitely small timestep, which we denote by ``dt``, from whence we see 
```math
P(dt) = \mathbb{I} + Q  dt
```
where ``\mathbb{I}`` is the identity matrix. Since the sum of the columns on the left hand side must be one and the sum of the columns of the identity matrix is one, we see that the only option is for the sum of the columns of ``Q`` to be zero for the identity to hold. The positivity of the off-diagonal terms of the transfer operator implies that the off-diagonal terms of the generator matrix, ``Q`` must be positive since the formula holds for every entry of the matrix and the identity matrix only modifies the diagonal terms. In fact we can say a bit more. Since probabilities are bounded above by one, we can say that the off-diagonal terms of the generator matrix must be less than or equal to zero.

## Using the generator to contruct Markov chains

We can now use the generator to construct Markov chains. The API in MarkovChainHammer.jl allows for two possibilites. The first is to construct the transfer operator from the generator and then directly use transfer operator as was done in the [Transfer Operators and Markov Chains](@ref sec:transfer_operators) section. 

For example, defining a generator ``Q`` and transfer operator ``P(\tau)`` as 

```@example generator2
using Random #hide
Random.seed!(1234) #hide
Q = [-1.0 2.0; 1.0 -2.0]
P(τ) = exp(Q * τ)
nothing #hide
```

We can define the ``\tau = 1``, transfer operator 

```@example generator2
P1 = P(1.0)
```

and then use the generate function as before, 
    
```@example generator2
using MarkovChainHammer.Trajectory: generate
markov_chain = generate(P1, 10)'
```

The second possibility is to use the generator directly. This is done by defining a generator ``Q``, choosing a step size ``dt``, and then using the ``generate`` function with the generator 

```@example generator2
n = 10; dt = 1.0;
markov_chain = generate(Q, n; dt = dt)'
```

The above is a chain of length 10, whose time intervals correspond to steps of size ``dt = 1.0``. The ``dt`` term is a keyword argument to the function and associated with a step size. Thus time "time" associated with the chain is, 
```@example generator2
time = dt * collect(0:n-1)'
```

We can also generate a markov chain sequence of length ``n=100`` with a step size ``dt = 0.1`` and then use the ``generate`` function to generate a chain of length ``n`` with a step size ``dt``. 

```@example generator2
n = 100; dt = 0.1;
markov_chain = generate(Q, n; dt = dt)'
```

where now the times associated with the process are 

```@example generator2
time = dt * collect(0:n-1)'
```

The flexibility of choosing the step size corresponding to a time, and still having the resulting process be [Markovian](https://en.wikipedia.org/wiki/Markov_property), is a powerful feature of the generator. Choosing an arbitrary ``dt`` but having the statistics remain the same over a physical timescale is what makes the resulting process a "continuous time markov process".

