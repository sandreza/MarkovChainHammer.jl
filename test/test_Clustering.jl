using MarkovChainHammer, Test, Revise, Random, LinearAlgebra

import MarkovChainHammer.Clustering: modularity_matrix, principal_vector, modularity
import MarkovChainHammer.Clustering: split_community, leicht_newman

@testset "Leicht Newman" begin
    # 3 communities
    A = [0 1 1; 1 0 1; 1 1 0]
    B = modularity_matrix(A)
    @test all(B .== [-1 -1 -1; -1 -1 -1; -1 -1 -1])
    s = principal_vector(B)
    @test all(s .== [-1, -1, 1])
    @test modularity(B) == -1
    community_cut1 = split_community(B, [1, 2, 3])
    @test all(community_cut1[1] .== [1])
    @test all(community_cut1[2] .== [2, 3])
    communities = leicht_newman(A)
    @test length(communities) == 3

    # 2 communities
    A = [1 1 0; 1 1 0; 0 0 1]
    B = modularity_matrix(A)
    @test all(B .≈ [0.4 -1.6 -2.4; -1.6 0.4 -2.4; -2.4 -2.4 2.4])
    s = principal_vector(B)
    @test all(s .== [-1, -1, 1])
    @test modularity(B) ≈ 9.6
    community_cut1 = split_community(B, [1, 2, 3])
    @test all(community_cut1[1] .== [3])
    @test all(community_cut1[2] .== [1, 2])
    communities = leicht_newman(A)
    @test length(communities) == 3

end