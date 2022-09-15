using OscoNet, Random, StableRNGs

@testset "Construct edge network" begin
    rng = StableRNG(123)
    n_permutations = 100
    data, _, _ = simulate_data(rng)
    adjacency_matrix, _, cost = bootstrap_hypothesis_test(rng, data, n_permutations)
 
    gene_names = ["gene_$i" for i = 1:20]
    edge_network = create_edge_network(adjacency_matrix, 1 ./ cost, gene_names)

    @test size(edge_network) == (47, 3)
    @test edge_network[45, 3] ≈ 0.03783900747679427
    @test edge_network[7, 3] ≈ 0.038445912047386455
end
