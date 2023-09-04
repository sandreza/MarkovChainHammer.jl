module MarkovChainHammer

using Reexport, PrecompileTools
@reexport using Statistics

include("TransitionMatrix/TransitionMatrix.jl")
include("Trajectory/Trajectory.jl")
include("BayesianMatrix/BayesianMatrix.jl")

include("Clustering/Clustering.jl")
include("Utils/Utils.jl")

@setup_workload begin
    using Random
    Random.seed!(1234)
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    markov_chain = rand([1, 2, 3], 100)
    @compile_workload begin
        QB = BayesianGenerator(markov_chain)
        QB = BayesianGenerator(markov_chain, QB.posterior)
        Qr = rand(QB)
        Qm = mean(QB)
        Q = generator(markov_chain)
        P = perron_frobenius(markov_chain)
        p = steady_state(Q)
        W = koopman_modes(Q)
        τs = decorrelation_times(Q)
        ht = holding_times(markov_chain)
        markov_chain_2 = generate(Q, 100)
    end
end

end # module MarkovChainHammer
