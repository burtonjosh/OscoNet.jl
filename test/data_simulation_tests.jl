using OscoNet

@testset "Simulate data" begin

    data, Ψ_g, ω = simulate_data()

    @test size(data) == (20, 1000)
    @test length(Ψ_g) == length(ω) == 20

    half_group = 80
    n_genes = 500
    n_cells = 500
    n_groups = 3

    data, Ψ_g, ω = simulate_data(;
        half_group = half_group,
        n_genes = n_genes,
        n_cells = n_cells,
        n_groups = n_groups,
    )

    @test size(data) == (500, 500)
    @test length(Ψ_g) == length(ω) == 500
end

@testset "True adjacency matrix" begin
    _, _, ω = simulate_data()
    true_matrix = true_adj_matrix(ω)

    @test size(true_matrix) == (20, 20)
    @test sum(true_matrix) == 90
end
