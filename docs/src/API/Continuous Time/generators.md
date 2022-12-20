# [Generators](@id sec:generators) 

In this section we discuss continuous time markov processes and methods for generating Markov chains from them. They are similar to a discrete time Markov chain, but with a notion of time built into their construction. The idea is to break up the time interval into a sequence of discrete time steps and then construct a discrete time Markov chain from the continuous time Markov process.

We start by building from the knowledge of the discrete time Markov chain and [transfer operators](@ref sec:transfer_operators). As mentioned before the columns of the transfer operator give information on where to go at a subsequent discrete step. But now we want to shift and develop a method that is consistent with the notion of continuous time. Thus the first thing we need to think about is, "How far do we go in time?" and then secondly "How do we make this consistent with the notion of a Markov process?". 

To motivate the use of a generator, we reexamine the transfer operator

```@example generator
ℳ = [0.6 0.3; 0.4 0.7]
ℳ
```

The transfer operator determines the probability of state j (here j = 1 or 2) to move to state i in one step. But suppose that we want to know the probability of moving from one state to another in *two* steps. Suppose that we are in state ```s_0`` at time 0. For example, if we were in state 2, then 
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
To get that 48% of our states ended up in state 1 and 52% of our states ended up in state 2. In total, if we want to just jump straight from our initial state ```s_0``` to state ```s_2```, we can just combine the two steps to get 
```@example generator
ℳ * ℳ * s₀
```

Thus we see if we want to move one step, we just multiply ```s_0``` by ```ℳ```, if we want to move two steps we multiply by, ```ℳ^2```. In general, if we want to move ```n```-steps, where ```n``` is a natural number, then we mulitply by ```ℳ^n```. If we want to move "0" steps, we can use the convection that ```ℳ^0``` is the identity matrix.

But now suppose that we don't want to compute at just integer steps, we want to compute what "happens in between", that is to say, raise to a non-integer power. One way to define this for matrices is to make use of the identity
```
\mathcal{M}^n = exp(n log(\mathcal{M}))
```
where the log and exp functions are the "matrix logarithm" and "matrix exponential". For example, the matrix logarithm of ```ℳ``` is
```@example generator
log(ℳ)
real.(log(ℳ)) #hide
``` 

and we can check 
```@example generator
exp(real.(log(ℳ))) - ℳ
real.(exp(real.(log(ℳ))) - ℳ) # hide
``` 

We can now check the formula
```@example generator
n = 2
exp(n * log(ℳ)) - ℳ^n
real.(exp(n * log(ℳ)) - ℳ^n) # hide
``` 

and finally, we can see what happens if we take non-integer steps
```@example generator
n = 1.5
ℳhalf = exp(n * log(ℳ)) 
``` 

and check to see taking two steps of size 1.5 is the same as taking one step of size 3,
```@example generator
ℳhalf - ℳ^3
``` 

This brings us to the point we were can finally introduce the generator. Presently, the generator is the matrix logarithm of the transfer operator 

```@example generator
Q = log(ℳ)
Q = real.(log(ℳ)) # hide
``` 

However, we have actually presented the relation backwards! The transfer operator for a time ```\tau``` in the future is the matrix exponential of the generator Q

```@example generator
P(τ) = exp(Q * τ)
``` 

In order for the formula to hold for all times ```\tau```, stricter requirements must be imposed on the generator. In order for ```P(\tau)``` to be a probability for all times ```\tau```, we require 
0. The sum of the columns of Q must be zero 
0. The off diagonal terms of Q must be greater than or equal to zero

These two requirements can be derived by taking ```\tau``` to be an infinitely small timestep, from whence we see 
```
P(dt) = I + Q * dt
```
where ```I``` is the identity matrix. Since the sum of the columns on the left hand side must be one and the sum of the columns of the identity matrix is one, we see that the only option is for the sum of the columns of ```Q``` to be zero for the identity to hold. The positivity of the off-diagonal terms of the transfer operator implies that the off-diagonal terms of the generator matrix, ```Q``` must be positive since the formula holds for every entry of the matrix and the identity matrix only modifies the diagonal terms.
