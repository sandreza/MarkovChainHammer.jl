using LinearAlgebra
function steady_state(Q)
    Λ, V = eigen(Q)
    return real(V[:, end] ./ sum(V[:,end]))
end

function entropy(p)
    N = length(p)
    entropy_value = 0.0
    for i in 1:N
        if p[i] > 100 * eps(eltype(p))
            entropy_value += p[i] .* log.(p[i])
        end
    end
    return -entropy_value / log(N)
end

function koopman_modes(Q)
    Λ, V = eigen(Q)
    V⁻¹ = inv(V)
    return V⁻¹
end

function decorrelation_times(Q)
    Λ = eigvals(Q)
    return 1 ./ real.(Λ)
end