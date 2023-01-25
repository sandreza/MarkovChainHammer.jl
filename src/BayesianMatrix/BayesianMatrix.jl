# Dependencies
# 

module BayesianMatrix

using Distributions, LinearAlgebra, Random, Statistics

# extend Base
import Base.rand
import LinearAlgebra.eigen, LinearAlgebra.eigvals, LinearAlgebra.size
import Statistics.mean, Statistics.var
import Base.show

# extend MarkovChainHammer 
import MarkovChainHammer.TransitionMatrix: holding_times, count_operator

# export functions
export BayesianGenerator, GeneratorParameterDistributions
export eigen_distribution, eigvals_distribution

# general abstractions
struct BayesianGenerator{PB,D,PA,PP}
    prior::PB
    data::D
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
