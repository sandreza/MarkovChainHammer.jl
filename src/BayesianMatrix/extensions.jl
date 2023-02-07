# Base.Random extensions
# extend Base
import Base.rand

function rand(Q::GeneratorParameterDistributions)
    rates = rand.(Q.rates)
    exit_probabilities = rand.(Q.exit_probabilities)
    return construct_generator(rates, exit_probabilities)
end

rand(Q::BayesianGenerator) = rand(Q.posterior)
rand(Q::GeneratorParameterDistributions, n::Int) = [rand(Q) for i in 1:n]
rand(Q::BayesianGenerator, n::Int) = [rand(Q) for i in 1:n]

# Base.Statistics 
import Statistics.mean, Statistics.var, Statistics.std
mean(Q::BayesianGenerator) = mean(Q.posterior)
mean(Q::GeneratorParameterDistributions) = construct_generator(mean.(Q.rates), mean.(Q.exit_probabilities))

function var(Q::GeneratorParameterDistributions)
    number_of_states = length(Q.rates)
    varQ = zeros(number_of_states, number_of_states)
    rates = Q.rates
    exit_probabilities = Q.exit_probabilities
    [varQ[i, i] = var(rates[i]) for i in eachindex(rates)]
    for i in 1:number_of_states
        # variance of product of two independent random variables
        σ₁² = var(rates[i])
        σ₂² = var(exit_probabilities[i])
        μ₁ = mean(rates[i])
        μ₂ = mean(exit_probabilities[i])
        varQ[[1:i-1..., i+1:number_of_states...], i] .= σ₁² * μ₂ .^ 2 + σ₂² * μ₁^2 + σ₁² * σ₂²
    end
    return varQ
end
var(Q::BayesianGenerator) = var(Q.posterior)

std(Q::BayesianGenerator) = std(Q.posterior)
std(Q::GeneratorParameterDistributions) = sqrt.(var(Q))

# Base.LinearAlgebra
import LinearAlgebra.eigen, LinearAlgebra.eigvals, LinearAlgebra.size

size(Q::BayesianGenerator) = (length(Q.prior.rates), length(Q.prior.rates))

eigen(Q::GeneratorParameterDistributions) = eigen(mean(Q))
eigen(Q::BayesianGenerator) = eigen(Q.posterior)

eigvals(Q::GeneratorParameterDistributions) = eigvals(mean(Q))
eigvals(Q::BayesianGenerator) = eigvals(Q.posterior)

eigen_distribution(Q::GeneratorParameterDistributions; samples=100) = eigen.(rand(Q, samples))
eigen_distribution(Q::BayesianGenerator; samples=100) = eigen_distribution(Q.posterior; samples=samples)

eigvals_distribution(Q::GeneratorParameterDistributions; samples=100) = eigvals.(rand(Q, samples))
eigvals_distribution(Q::BayesianGenerator; samples=100) = eigvals_distribution(Q.posterior; samples=samples)

# Base.copy
import Base.copy
copy(Q::GeneratorParameterDistributions) = GeneratorParameterDistributions(copy(Q.rates), copy(Q.exit_probabilities))
copy(Q::GeneratorPredictiveDistributions) = GeneratorPredictiveDistributions(copy(Q.holding_times), copy(Q.exit_counts))
copy(Q::BayesianGenerator) = BayesianGenerator(copy(Q.prior), copy(Q.data), copy(Q.posterior), copy(Q.predictive))

# Base.show
import Base.show

function Base.show(io::IO, Q::GeneratorParameterDistributions)
    color1 = :white
    color2 = :magenta
    color3 = :cyan

    printstyled(io, "GeneratorParameterDistributions \n", color=color1)
    printstyled(io, "rates  \n", color=color3)
    show(io, MIME("text/plain"), Q.rates)
    println(" ")
    printstyled(io, "exit_probabilities  \n", color=color2)
    show(io, MIME("text/plain"), Q.exit_probabilities)
    return nothing
end

function Base.show(io::IO, Q::BayesianGenerator)
    color1 = :white
    color2 = :magenta
    color3 = :cyan
    printstyled(io, "BayesianGenerator \n", color=color1)
    printstyled(io, "prior variance  \n", color=color2)
    show(io, MIME("text/plain"), var(Q.prior))
    println(" ")
    printstyled(io, "prior mean  \n", color=color2)
    show(io, MIME("text/plain"), mean(Q.prior))
    println(" ")
    printstyled(io, "posterior variance  \n", color=color3)
    show(io, MIME("text/plain"), var(Q.posterior))
    println(" ")
    printstyled(io, "posterior mean  \n", color=color3)
    show(io, MIME("text/plain"), mean(Q.posterior))

    return nothing
end