using OscoNet

@testset "Construct edge network" begin
    Random.seed!(1234)
    n_permutations = 100
    data, _, _ = simulate_data()
    adjacency_matrix, _, cost = bootstrap_hypothesis_test(data, n_permutations)
    gene_names = ["gene_$i" for i = 1:20]

    edge_network = create_edge_network(adjacency_matrix, 1 ./ cost, gene_names)

    @test size(edge_network) == (45, 3)
    @test edge_network[45, 3] ≈ 0.0271483972
    @test edge_network[7, 3] ≈ 0.3718996217
end
