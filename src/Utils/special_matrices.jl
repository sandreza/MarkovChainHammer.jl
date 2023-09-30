export ornstein_uhlenbeck_generator, central_advection_periodic, discrete_laplacian_periodic

function ornstein_uhlenbeck_generator_raw(n)
    Mⱼₖ = zeros(n + 1, n + 1)
    δ(j, k) = (j == k) ? 1 : 0

    for j in 0:n, k in 0:n
        jj = j + 1
        kk = k + 1
        Mⱼₖ[jj, kk] =
            (-n * δ(j, k) + k * δ(j + 1, k) + (n - k) * δ(j - 1, k)) / 2
    end
    return Mⱼₖ
end

ornstein_uhlenbeck_generator(n) = ornstein_uhlenbeck_generator_raw(n-1)


function central_advection_periodic(N; Δx=1.0, u=1.0)
    A = zeros(N, N)
    for i in 1:N
        A[(i-1)%N+1, (i-1*0)%N+1] = -1
        A[(i+1)%N+1, (i-1*0)%N+1] = 1
    end
    return u * A / 2Δx
end

function discrete_laplacian_periodic(N; Δx=1.0, κ=1.0)
    Δ = zeros(N, N)
    for i in 1:N
        Δ[i, i] = -2
        Δ[(i-1)%N+1, (i-1*0)%N+1] = 1
        Δ[(i+1)%N+1, (i-1*0)%N+1] = 1
    end
    return κ * Δ / (Δx^2)
end
