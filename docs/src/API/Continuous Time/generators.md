# [Generators](@id sec:generators) 

In this section we discuss continuous time markov processes and methods for generating Markov chains from them. They are similar to a discrete time Markov chain, but with a notion of time built into their construction. The idea is to break up the time interval into a sequence of discrete time steps and then construct a discrete time Markov chain from the continuous time Markov process.

We start by building from the knowledge of the discrete time Markov chain and [transfer operators](@ref sec:transfer_operators). As mentioned before the columns of the transfer operator give information on where to go at a subsequent discrete step. But now we want to shift and develop a method that is consistent with the notion of continuous time. Thus the first thing we need to think about is, "How far do we go in time?" and then secondly "How do we make this consistent with the notion of a Markov process?". 

