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
function GeneratorParameterDistributions(number_of_states::Int;  α=ones(number_of_states), β=ones(number_of_states), αs=ones(number_of_states - 1, number_of_states))
    if length(α) == 1 && length(β) == 1 
        α = fill(α, number_of_states)
        β = fill(β, number_of_states)
        αs = repeat(αs, 1, number_of_states)
    end
    @assert size(αs)[1] == number_of_states - 1
    @assert size(αs)[2] == size(α)[1] == size(β)[1] == number_of_states
    @assert all(αs .> 0)
    @assert all(α .> 0)
    @assert all(β .> 0)
    α = α
    β = β
    θ = 1 ./ β
    rates = [Gamma(α[i], θ[i]) for i in 1:number_of_states]
    exit_probabilities = [Dirichlet(αs[:,i]) for i in 1:number_of_states]
    return GeneratorParameterDistributions(rates, exit_probabilities)
end

function GeneratorParameterDistributions(parameters::NamedTuple)
    @assert haskey(parameters, :α)
    @assert haskey(parameters, :β)
    @assert haskey(parameters, :αs)
    @assert length(parameters.α) == length(parameters.β) == length(parameters.αs[1,:]) == (length(parameters.αs[:,1])+1)
    return GeneratorParameterDistributions(length(parameters.α), α=parameters.α, β=parameters.β, αs=parameters.αs)
end

function GeneratorPredictiveDistributions(number_of_states::Int; μ=ones(number_of_states), σ=ones(number_of_states), ξ=ones(number_of_states), n=ones(number_of_states), αs=ones(number_of_states - 1, number_of_states))
    if length(μ) == 1 && length(σ) == 1
        μ = fill(α, number_of_states)
        σ = fill(β, number_of_states)
        ξ = fill(ξ, number_of_states)
        n = fill(n, number_of_states)
        αs = repeat(αs, 1, number_of_states)
    end
    @assert size(αs)[1] == number_of_states - 1
    @assert size(αs)[2] == size(μ)[1] == size(ξ)[1] == size(n)[1] == number_of_states
    @assert all(αs .> 0)
    @assert all(n .> 0)
    @assert all(σ .> 0)

    holding_times = [GeneralizedPareto(μ[i], σ[i], ξ[i]) for i in 1:number_of_states]
    exit_counts = [DirichletMultinomial(n[i], αs[:,i]) for i in 1:number_of_states]
    return GeneratorPredictiveDistributions(holding_times, exit_counts)
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
    p_data = count_operator(data, number_of_states)
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
            @warn "no data for state $i, falling back on prior"
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
        αs_new = αs + p_data[[1:i-1..., i+1:number_of_states...], i]
        posterior_exit_probabilities[i] = Dirichlet(αs_new)
        # predictive
        if number_of_exits > 0
            n = length(ht_data[i])
        else
            n = 1
        end
        predictive_exit_counts[i] = DirichletMultinomial(n, αs_new)
    end
    posterior = GeneratorParameterDistributions(posterior_rates, posterior_exit_probabilities)
    predictive = GeneratorPredictiveDistributions(predictive_holding_times, predictive_exit_counts)
    return BayesianGenerator(prior, posterior, predictive)
end

BayesianGenerator(prior::GeneratorParameterDistributions) = BayesianGenerator(prior, prior, nothing)
BayesianGenerator(data; dt = 1) = BayesianGenerator(data,  uninformative_prior(maximum(data)) ; dt=dt)

"""
    BayesianGenerator(data::Vector{Vector{Int64}}, prior::GeneratorParameterDistributions; dt=1)

    The ensemble version of the BayesianGenerator constructor. Here the data is a vector of vectors of integers, where each vector of integers is a single ensemble member trajectory.

    # Arguments
    - `data::Vector{Vector{Int64}}`: The data to construct the BayesianGenerator object from.
    - `prior::GeneratorParameterDistributions`: The prior distribution for the BayesianGenerator object.

    # Keyword Arguments
    - `dt::Number=1`: The time step between each data point.

    # Returns
    - `BayesianGenerator`: A BayesianGenerator object constructed from the data and the prior distribution. Contains the posterior distributions for the rates and exit probabilities, as well as the predictive distributions for the holding times and exit counts.
"""
function BayesianGenerator(data::Vector{Vector{Int64}}, prior::GeneratorParameterDistributions; dt=1)
    new_prior = prior
    Q = BayesianGenerator(data[1], new_prior; dt=dt)
    for i in eachindex(data)
        @assert all(data[i] .> 0)
        Q = BayesianGenerator(data[i], new_prior; dt=dt)
        new_prior = Q.posterior
    end
    return Q
end
BayesianGenerator(data::Vector{Vector{Int64}}; dt=1) = BayesianGenerator(data, uninformative_prior(maximum(reduce(vcat, data))); dt=dt)

"""
    BayesianGenerator(parameters::NamedTuple)

Construct a BayesianGenerator object from a NamedTuple of parameters.

# Arguments
- `parameters::NamedTuple`: A NamedTuple containing the parameters for the BayesianGenerator object. Must contain the fields `prior`, `posterior` and `predictive`, each of which must be a NamedTuple containing the parameters for the respective distributions.

# Returns
- `BayesianGenerator`: A BayesianGenerator object constructed from the parameters. Contains the posterior distributions for the rates and exit probabilities, as well as the predictive distributions for the holding times and exit counts.

"""
function BayesianGenerator(parameters::NamedTuple) 
    @assert haskey(parameters, :prior)
    @assert haskey(parameters, :posterior)
    @assert haskey(parameters, :predictive)

    number_of_states = length(parameters.prior.α)
    prior = GeneratorParameterDistributions(number_of_states; α=parameters.prior.α, β=parameters.prior.β, αs=parameters.prior.αs)
    posterior = GeneratorParameterDistributions(number_of_states; α=parameters.posterior.α, β=parameters.posterior.β, αs=parameters.posterior.αs)
    predictive = GeneratorPredictiveDistributions(number_of_states; μ=parameters.predictive.μ, σ=parameters.predictive.σ, ξ=parameters.predictive.ξ, n=parameters.predictive.n, αs=parameters.predictive.αs)
    return BayesianGenerator(prior, posterior, predictive)
end

"""
    uninformative_prior(number_of_states; scale=eps(1e2))

Construct an uninformative prior distribution for a BayesianGenerator object.

# Arguments
- `number_of_states::Int`: The number of states in the BayesianGenerator object.

# Keyword Arguments

- `scale::Number=eps(1e2)`: The scale parameter for the Gamma and Dirichlet distributions.

# Returns
- `GeneratorParameterDistributions`: An uninformative prior distribution for a BayesianGenerator object.

"""
uninformative_prior(number_of_states; scale=eps(1e2)) = GeneratorParameterDistributions(number_of_states::Int; α=scale, β=scale, αs=ones(number_of_states - 1) * scale)

"""
    unpack(Q::BayesianGenerator)

Unpack a BayesianGenerator object into its parameters.

# Arguments
- `Q::BayesianGenerator`: The BayesianGenerator object to unpack.

# Returns
- `Tuple{NamedTuple{(:prior, :posterior, :predictive),Tuple{NamedTuple{(:α, :β, :αs),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}},NamedTuple{(:α, :β, :αs),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}},NamedTuple{(:μ, :σ, :ξ, :n, :αs),Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Vector{Int64},Vector{Vector{Float64}}}}}}}`: A tuple containing the parameters for the prior, posterior and predictive distributions.

"""
function unpack(Q::BayesianGenerator)
    prior = params(Q.prior)
    posterior = params(Q.posterior)
    predictive = params(Q.predictive)
    return (; prior, posterior, predictive)
end
