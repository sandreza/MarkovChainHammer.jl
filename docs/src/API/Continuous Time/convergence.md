# [Convergence of Generators](@id sec:generator_convergence)

Convergence of a generator from data is more subtle than the convergence of a transfer operator from data due to the presence of time scales in the former.  For a fixed timestep size ``\Delta t`` and infinite data we have convergence to the associated transfer operator. Thus we take the limit of ``\Delta t \rightarrow 0`` afterwards to get convergence to the generator. 

Given the exact generator

```@example generator_convergence 
Q_exact = [-1.0 2.0; 1.0 -2.0]
Q_exact
```

We generate a Markov chains of increasing size at a fixed timestep ``dt`` and compute the empirical transition matrix. We then verify that the empirical transition matrix converges to the exact transition operator. Let us first generate a Markov chain of size 100.

```@example generator_convergence
using MarkovChainHammer.Trajectory: generate
using Random #hide
Random.seed!(1234) #hide
dt = 1.0
steps = 100
markov_chain = generate(Q_exact, steps; dt=dt); 
nothing; # hide
```

We now load the empirical operator constructors for both the generator and the transfer operator in order to seperately check for convergence. 

First let us construct the transfer operator
```@example generator_convergence
using MarkovChainHammer.TransitionMatrix: perron_frobenius
ℳ_empirical = perron_frobenius(markov_chain)
```

We see that the above matrix is close to the exact transfer operator
```@example generator_convergence
ℳ_exact = exp(Q_exact)
```

Now let us construct the generator
```@example generator_convergence
using MarkovChainHammer.TransitionMatrix: generator
Q_empirical = generator(markov_chain; dt = dt)
```

We see that the above matrix is quite far from the exact generator, which leads to the question "what happened?". The resolution to this dilemma is to observe that timescales of the markov chain were too large as compared to the expected holding times of the exact generator, which are 1.0 for state 1 and ``1/2`` for state 2. If the chain is instead generated with smaller timesteps 

```@example generator_convergence
dt = 0.1
markov_chain = generate(Q_exact, steps; dt=dt); 
Q_empirical = generator(markov_chain; dt = dt)
```

We can see that the estimate of the generator has improved. Note that we can also calculate the generator by taking the matrix logarithm and dividing by the timescale to get a similar answer

```@example generator_convergence 
log(perron_frobenius(markov_chain)) / dt
```


Let us now automate the process of checking convergence as function of fixed timestep size. We generate a Markov chain of increasing size and compute the empirical transition matrix. We then verify that the empirical transition matrix converges to the exact transition operator. We use the function `norm` from the `LinearAlgebra` package to compute the norm of the difference between the empirical and exact transition matrices. We use the difference between the matrices divided by the norm of the exact matrix to get a relative error. At small timesteps the transfer operator converges to the identity matrix, hence we divide the error by a factor of ``dt`` in order to make a fair comparison to the relative error of the generator matrix. The takeaway message is that the error of the generator is bounded by both the timestep size and the available number of independent samples in the data. Here is the code to verify this statement

```@example generator_convergence
using LinearAlgebra: norm
println("Transfer Operator Convergence")
for dt in [1.0, 0.1, 0.01, 0.001]
    ℳ_exact_local = exp(Q_exact * dt)
    println("For a timestep of $dt")
    for i in 2:2:6
        n = 10^i
        markov_chain_local = generate(ℳ_exact_local, n)
        number_of_states = 2
        ℳ_empirical_local = perron_frobenius(markov_chain_local, number_of_states)
        empirical_error = norm(ℳ_empirical_local - ℳ_exact_local) / norm(ℳ_exact_local) / dt
        println("A chain of size 10^$i yields a relative empirical error of $(empirical_error)")
    end
    println("---------------------------")
end

println("---------------------------")
println("Generator Convergence")
for dt in [1.0, 0.1, 0.01, 0.001]
    ℳ_exact_local = exp(Q_exact * dt)
    println("For a timestep of $dt")
    for i in 2:2:6
        n = 10^i
        markov_chain_local = generate(ℳ_exact_local, n)
        number_of_states = 2
        Q_empirical_local = generator(markov_chain_local, number_of_states; dt = dt)
        empirical_error = norm(Q_empirical_local - Q_exact) / norm(Q_exact)
        println("A chain of size 10^$i yields a relative empirical error of $(empirical_error)")
    end
    println("---------------------------")
end
```

In particular, for the convergence of the generator, we see that the error is roughly on the order of the timestep size, except for the last case. A heuristic explanation for this behavior is as follows. If we think of each entry of the matrix as a random variable then the convergence to the expected value (the exact entry of the matrix) is inversely proportional to the square root of the number of independent samples. As one decreases the timestep size, the number of independent samples decrease since the number of timesteps is fixed, that is to say we are looking at the simulation for an increasingly shorter amount of physical time. Roughly speaking one would expect the error to be 

```math
error \propto \max \left(dt, \frac{1}{\sqrt{S}}  \right)
```

where ``S`` is the number of independent samples. In the last case, since the largest decorrelation time of the markov chain is on the order of 1 time unit, we expect the number of independent samples in a physical time interval of ``10^6  dt = 10^3`` to be roughly on the order of ``10^3 / 5``, where 5 time units is taken as a decorrelation threshold for the markov chain. We see that the error is 
```@example generator_convergence
approximate_error = 1 / sqrt(10^3 / 5)
```
which is roughly the same order of magnitude as the final error.

