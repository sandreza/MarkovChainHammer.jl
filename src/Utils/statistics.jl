

function autocovariance(x; timesteps=length(x))
    μ = mean(x)
    autocor = zeros(timesteps)
    for i in ProgressBar(1:timesteps)
        autocor[i] = mean(x[i:end] .* x[1:end-i+1]) - μ^2
    end
    return autocor
end

function autocovariance(g⃗, Q::Eigen, timelist)
    autocov = zeros(length(timelist))
    # Q  = V Λ V⁻¹
    Λ, V = Q
    p = real.(V[:, end] ./ sum(V[:, end]))
    v1 = V \ (p .* g⃗)
    w1 = g⃗' * V
    μ = sum(p .* g⃗)
    for i in eachindex(timelist)
        autocov[i] = real(w1 * (exp.(Λ .* timelist[i]) .* v1)) - μ^2
    end
    return autocov
end

function autocovariance(observable, Ps::Vector{Matrix{Float64}}, steps)
    autocor = zeros(steps + 1)
    p = steady_state(Ps[1])
    μ² = sum(observable .* p)^2
    autocor[1] = observable' * (observable .* p) - μ²
    for i in 1:steps
        autocor[i+1] = observable' * Ps[i] * (observable .* p) - μ²
    end
    return autocor
end

autocovariance(g⃗, Q, timelist) = autocovariance(g⃗, eigen(Q), timelist)
