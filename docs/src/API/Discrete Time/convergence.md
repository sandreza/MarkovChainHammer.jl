# [Convergence of Transfer Operators](@id sec:to_convergence)

In this section we show how to verify convergence of transfer operators. We use the following example: 

Given the exact transfer operator

```@example transfer_operator_convergence 
ℳ_exact = [0.6 0.3; 0.4 0.7]
ℳ_exact
```

We generate a Markov chains of increasing size and compute the empirical transition matrix. We then verify that the empirical transition matrix converges to the exact transition operator. Let us first generate a Markov chain of size 100.

```@example transfer_operator_convergence
using MarkovChainHammer.Trajectory: generate
using Random #hide
Random.seed!(1234) #hide
markov_chain = generate(ℳ_exact, 100); 
nothing; # hide
```

and then compute the empirical operator
```@example transfer_operator_convergence
using MarkovChainHammer.TransitionMatrix: perron_frobenius
ℳ_empirical = perron_frobenius(markov_chain)
ℳ_empirical
```

We see that the empirical transition matrix is close to the exact transition matrix. Let us now generate a Markov chain of size 10000 and compute the empirical transition matrix.

```@example transfer_operator_convergence
markov_chain = generate(ℳ_exact, 10000);
ℳ_empirical = perron_frobenius(markov_chain)
ℳ_empirical
```

We see that the empirical transition matrix is closer to the exact transition matrix. 

Let us now automate the process of checking convergence. We generate a Markov chain of increasing size and compute the empirical transition matrix. We then verify that the empirical transition matrix converges to the exact transition operator. We use the function `norm` from the `LinearAlgebra` package to compute the norm of the difference between the empirical and exact transition matrices.

```@example transfer_operator_convergence
using LinearAlgebra: norm
for i in 1:6
    n = 10^i
    markov_chain_local = generate(ℳ_exact, n)
    ℳ_empirical_local = perron_frobenius(markov_chain_local)
    empirical_error = norm(ℳ_empirical_local - ℳ_exact)
    println("A chain of size 10^$i yields an empirical error of $(empirical_error)")
end
```

We see that the error in the empirical estimate decreases as the chain size increases.

Here we used knowledge of the exact operator to show that the operator converges. If we only have the data, as is often the case in practice, then we cannot know if the empirical operator is sufficiently well-estimated; however, sometimes we know particular abstract properties of the Markov chain which imply theoretical gaurantees on the convergence of the empirical operator in the limit of infinite data. On the other end, given a finite time-series, we can view the operator as a compression algorithm for the stochastic process in which case we only look for data-consistency.  