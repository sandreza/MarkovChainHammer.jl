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


Let us now automate the process of checking convergence as function of fixed timestep size. We generate a Markov chain of increasing size and compute the empirical transition matrix. We then verify that the empirical transition matrix converges to the exact transition operator. We use the function `norm` from the `LinearAlgebra` package to compute the norm of the difference between the empirical and exact transition matrices.

```@example generator_convergence
using LinearAlgebra: norm
for dt in [1.0, 0.1, 0.01]
    ℳ_exact = exp(Q_exact * dt)
    for i in 2:2:6
        n = 10^i
        markov_chain_local = generate(ℳ_exact, n)
        ℳ_empirical_local = perron_frobenius(markov_chain_local)
        empirical_error = norm(ℳ_empirical_local - ℳ_exact) / norm(ℳ_exact)
        println("A chain of size 10^$i with timestep $dt yields a relative empirical error of $(empirical_error) for the transfer operator")
    end
    println("---------------------------")
end

for dt in [1.0, 0.1, 0.01]
    ℳ_exact = exp(Q_exact * dt)
    for i in 2:2:6
        n = 10^i
        markov_chain_local = generate(ℳ_exact, n)
        Q_empirical_local = generator(markov_chain_local; dt)
        empirical_error = norm(Q_empirical_local - Q_exact) / norm(Q_exact)
        println("A chain of size 10^$i with timestep $dt yields a relative empirical error of $(empirical_error) for the generator")
    end
    println("---------------------------")
end
```

The takeaway message in the present context is that the error of the generator is bounded by both the timestep size and the available number of independent samples in the data. 
