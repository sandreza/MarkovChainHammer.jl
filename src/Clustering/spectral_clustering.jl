using LinearAlgebra

# credit to Peter Schmid and Ludovico Giorgini

function modularity_matrix(A; ϵ=eps(100.0))
    abs_A = abs.(A)
    maxA = maximum(abs_A)
    Ã = (abs_A .>= ϵ * maxA) + I
    m = sum(A)
    Ki = sum(Ã, dims=1) # in-degree
    Ko = sum(Ã, dims=2) # out-degree
    b = Ã - (Ko * Ki) / m
    B = Symmetric(b + b')
    return B
end

function giorgini_matrix(A)
    N = size(A)[1]
    b = zeros(N, N)
    for i = 1:N, j = 1:N
        b[i, j] = (A[i, j] - (1 - A[i, i]) * (sum(A[j, :])) / N) / N
    end
    B = Symmetric(b + b')
    return B
end

function generator_fine_matrix(Q)
    Λ = eigenvalues(Q)
    P = exp(Q / real(Λ[1]))
    return giorgini_matrix(P)
end

function generator_coarse_matrix(Q)
    Λ = eigenvalues(Q)
    P = exp(Q / real(Λ[end-1]))
    return giorgini_matrix(P)
end

function principal_vector(B::Symmetric)
    s = ones(Int, size(B)[1])
    Λ, V = eigen(B)
    v₁ = V[:, argmax(Λ)]
    s[v₁.<=0] .= -1
    return s
end

function modularity(B, s)
    return s' * (B * s)
end

modularity(B::Symmetric) = modularity(B, principal_vector(B))

function split_community(B, indices; ϵ=eps(1e6))
    Bg = B[indices, :][:, indices]
    Bg = Bg - Diagonal(sum(Bg, dims=2)[:])
    Bg = Symmetric(Bg + Bg')
    s = principal_vector(Bg)
    q = modularity(Bg, s)

    if (q > ϵ)
        ind1 = [i for (j, i) in enumerate(indices) if s[j] == 1]
        ind2 = [i for (j, i) in enumerate(indices) if s[j] == -1]
        return ind1, ind2
    end

    return [], []

end

"""
`leicht_newman(A; modularity_type = :giorgini)`

# Description 
    Compute the communities of a graph using the Leicht Newman algorithm.

# Arguments
- `A::AbstractArray`: Adjacency matrix of the graph.

# Keyword Arguments
- `modularity_type::Symbol`: Type of modularity matrix to use. Can be `:giorgini`, `:modularity`, `:generator_fine`, or `:generator_coarse`. Defaults to `:modularity`.

# Returns
- `AbstractArray`: Array of communities.
"""
function leicht_newman(A; modularity_type = :modularity)
    if modularity_type == :giorgini
        B = giorgini_matrix(A)
    elseif modularity_type == :modularity
        B = modularity_matrix(A)
    elseif modularity_type == :generator_fine
        B = generator_fine_matrix(A)
    elseif modularity_type == :generator_coarse
        B = generator_coarse_matrix(A)
    else
        error("modularity_type must be :giorgini , :modularity, :generator_fine, or :generator_coarse")
    end
    n = size(A)[1]
    W, F = [collect(1:n)], []

    while (length(W) > 0)
        w = popfirst!(W)
        ind1, ind2 = split_community(B, w)
        if (length(ind1) > 0) & (length(ind2) > 0)
            W = [ind1, ind2, W...]
        else
            push!(F, w)
        end
    end

    return F
end