# utility function 
function construct_generator(rates, exit_probabilities)
    number_of_states = length(rates)
    generator = zeros(number_of_states, number_of_states)
    generator -= I
    for i in 1:number_of_states
        generator[[1:i-1..., i+1:number_of_states...], i] .= exit_probabilities[i]
        generator[:, i] *= rates[i]
    end
    return generator
end

"""
    GeneratorParameterDistributions(number_of_states::Int; α=1, β=1, αs=ones(number_of_states - 1))

Construct a GeneratorParameterDistributions object with Gamma(α, 1/β) prior distributions for the rates and Dirichlet(αs) prior distributions for the exit probabilities. Each state has the same probability distribution
This is useful for construction prior distributions quickly. The underlying assumption for default rates is that the units of time are 1, and thus the rates are all by default 1. Futhermore the mean probability of each state is  1/(number_of_states-1)
    
The reason why the αs is of length(number_of_states - 1) is because they are exit probabilities, and hence the probability of returning to the same state is 0.
The index ordering for the exit probability of state i, is [1:i-1..., i+1:number_of_states...]. 
For example, if number_of_states = 3, then the index ordering for state 1 is indices = [2, 3], i.e., indices[1] = 2 and indices[2] = 3, meaning, the first index corresponds to exiting through state 2 and the second index corresponds to exiting through state 3.
The index ordering for state 2 is indices = [1, 3], i.e., indices[1] = 1 and indices[2] = 3, meaning, the first index corresponds to exiting through state 1 and the second index corresponds to exiting through state 3.
The index ordering for state 3 is indices = [1, 2], i.e., indices[1] = 1 and indices[2] = 2, meaning, the first index corresponds to exiting through state 1 and the second index corresponds to exiting through state 2.

# Arguments
- `number_of_states::Int`: The number of states in the Markov chain.

# Keyword Arguments
- `α::Number=2`: The shape parameter for the Gamma prior distribution for the rates.
- `β::Number=2`: The rate parameter for the Gamma prior distribution for the rates.
- `αs::Vector{Number}=ones(number_of_states - 1)`: The concentration parameters for the Dirichlet prior distribution for the exit probabilities.

# Returns

- `GeneratorParameterDistributions`: A GeneratorParameterDistributions object with Gamma(α, 1/β) prior distributions for the rates and Dirichlet(αs) prior distributions for the exit probabilities.

# Note on the Gamma prior distribution (from https://en.wikipedia.org/wiki/Gamma_distribution)

The Gamma distribution is a conjugate prior distribution for the exponential family distribution of the rates.
The Gamma distribution is parameterized by the shape parameter α and the rate parameter β. The mean of the Gamma distribution is α/β and the variance is α/β^2.

# Note on the Dirichlet prior distribution (from https://en.wikipedia.org/wiki/Dirichlet_distribution)

The Dirichlet prior distribution is a conjugate prior distribution for the multinomial distribution of the exit probabilities.
The Dirichlet prior distribution is parameterized by the concentration parameters αs. 
The mean of the Dirichlet distribution is E[X⃗] = αs/sum(αs) 
The covariance is CoVar[X⃗ ⊗ X⃗] = Diagonal(α̃) - α̃ ⊗ α̃ / (α₀ + 1) where α̃ = αs / α₀  and α₀ = sum(αs).
The variance Var(Xᵢ) = α̃ᵢ(1-α̃ᵢ) / (α₀ + 1) where α₀ = sum(αs).

"""
function GeneratorParameterDistributions(number_of_states::Int; α=1, β=1, αs=ones(number_of_states - 1))
    @assert length(αs) == number_of_states - 1
    @assert all(αs .> 0)
    @assert α > 0
    @assert β > 0
    α = α
    β = β
    θ = 1 / β
    rates = [Gamma(α, θ) for i in 1:number_of_states]
    exit_probabilities = [Dirichlet(αs) for i in 1:number_of_states]
    return GeneratorParameterDistributions(rates, exit_probabilities)
end

"""
    BayesianGenerator(data, prior::GeneratorParameterDistributions; dt=1)

Construct a BayesianGenerator object from data and a prior distribution.

# Arguments
- `data::Vector{Int}`: The data to construct the BayesianGenerator object from.
- `prior::GeneratorParameterDistributions`: The prior distribution for the BayesianGenerator object.

# Keyword Arguments
- `dt::Number=1`: The time step between each data point.

# Returns
- `BayesianGenerator`: A BayesianGenerator object constructed from the data and the prior distribution. Contains the posterior distributions for the rates and exit probabilities, as well as the predictive distributions for the holding times and exit counts.

"""
function BayesianGenerator(data, prior::GeneratorParameterDistributions; dt=1)
    number_of_states = length(prior.rates)
    ht_data = holding_times(data, number_of_states; dt=dt)
    p_data = Int.(count_operator(data, number_of_states))
    p_data = p_data - Diagonal(diag(p_data))

    posterior_rates = Vector{Gamma{Float64}}(undef, number_of_states)
    posterior_exit_probabilities = Vector{Dirichlet{Float64,Vector{Float64},Float64}}(undef, number_of_states)
    predictive_holding_times = Vector{GeneralizedPareto{Float64}}(undef, number_of_states)
    predictive_exit_counts = Vector{DirichletMultinomial{Float64}}(undef, number_of_states)

    for i in 1:number_of_states
        number_of_exits = length(ht_data[i])
        α, θ = params(prior.rates[i])
        # first handle the rates
        # posterior
        β = 1 / θ
        α_new = α + number_of_exits
        if number_of_exits > 0
            β_new = β + sum(ht_data[i])
        else
            @warn "no data for state $i, falling back on prior for posterior distribution"
            β_new = β
        end
        θ_new = 1 / β_new
        posterior_rates[i] = Gamma(α_new, θ_new)
        # predictive posterior 
        μ = 0 # lomax
        σ = β_new / α_new
        ξ = 1 / α_new
        predictive_holding_times[i] = GeneralizedPareto(μ, σ, ξ)

        # next the exit probabilities
        # posterior
        αs = params(prior.exit_probabilities[i])[1]
        αs_new = αs + p_data[i, [1:i-1..., i+1:number_of_states...]]
        posterior_exit_probabilities[i] = Dirichlet(αs_new)
        # predictive
        if number_of_exits > 0
            n = length(ht_data[i])
        else
            @warn "no data for state $i, falling back on DirichletMultinomial with one observation for predictive distribution"
            n = 1
        end
        predictive_exit_counts[i] = DirichletMultinomial(n, αs_new)
    end
    posterior = GeneratorParameterDistributions(posterior_rates, posterior_exit_probabilities)
    predictive = GeneratorPredictiveDistributions(predictive_holding_times, predictive_exit_counts)
    return BayesianGenerator(prior, data, posterior, predictive)
end

BayesianGenerator(prior::GeneratorParameterDistributions) = BayesianGenerator(prior, nothing, prior, nothing)
BayesianGenerator(data; dt = 1) = BayesianGenerator(data,  uninformative_prior(maximum(data)) ; dt=dt)

uninformative_prior(number_of_states) = GeneratorParameterDistributions(number_of_states::Int; α=eps(10.0), β=eps(10.0), αs=ones(number_of_states - 1) * eps(10.0))