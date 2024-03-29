@reexport module BayesianMatrix

using Distributions, LinearAlgebra, Random, Statistics

# load MarkovChainHammer 
import MarkovChainHammer.TransitionMatrix: holding_times, count_operator

# export functions
export BayesianGenerator, GeneratorParameterDistributions
export eigen_distribution, eigvals_distribution
export uninformative_prior 
export params, unpack

# general abstractions
struct BayesianGenerator{PB,PA,PP}
    prior::PB
    posterior::PA
    predictive::PP
end

struct GeneratorParameterDistributions{H,P}
    rates::H
    exit_probabilities::P
end

struct GeneratorPredictiveDistributions{H,P}
    holding_times::H
    exit_counts::P
end

include("constructors.jl")

include("extensions.jl")

end # module BayesianMatrix
