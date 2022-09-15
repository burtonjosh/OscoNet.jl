using OscoNet, Random, LinearAlgebra, StableRNGs

@testset "Hypothesis test" begin
    rng = StableRNG(123)
    n_permutations = 100

    data, _, _ = simulate_data(rng)

    adjacency_matrix, qvalues, cost = bootstrap_hypothesis_test(rng, data, n_permutations)

    @test size(adjacency_matrix) == (20, 20)
    @test size(qvalues) == (20, 20)
    @test size(cost) == (20, 20)

    @test sum(adjacency_matrix) == 94
end

@testset "Ψ for each gene pair" begin
    rng = StableRNG(123)
    data, _, _ = simulate_data(rng)

    Ψ, cost = OscoNet.find_best_psi_for_each_gene_pair(data)

    @test sum(OscoNet.vec_triu_loop(Ψ)) ≈ 302.04322717500526
    @test sum(OscoNet.vec_triu_loop(cost)) ≈ 441848.65225983947
end

@testset "Permuted cost" begin
    n_permutations = 100
    data, _, _ = simulate_data()
    cost_permuted = OscoNet.get_permuted_cost(data, n_permutations)

    @test size(cost_permuted) == (20, 20, 100)
end

@testset "Matrix manipulation" begin
    A = [1 2 3; 4 5 6; 7 8 9]
    @test OscoNet.vec_triu_loop(A) == [2, 3, 6]

    # construct vector of matrix indices
    upper_triangle_indices = [20 * (j - 1) + i for i = 2:20 for j = 1:(i-1)]
    A = reshape(1:400, 20, 20)'

    @test OscoNet.vec_triu_loop(A) == upper_triangle_indices

    sym_A = OscoNet.construct_symmetric_matrix_from_upper_triangular_vector(
        upper_triangle_indices,
    )

    @test issymmetric(sym_A)
    @test size(sym_A) == (20, 20)
end

@testset "Triangular numbers" begin
    @test OscoNet.istriangular(45)
    @test OscoNet.istriangular(45.0)
    @test OscoNet.istriangular(0)
    @test OscoNet.istriangular(0.0)
    @test OscoNet.istriangular(780)
    @test OscoNet.istriangular(780.0)
    @test OscoNet.istriangular(2.0) == false

    @test OscoNet.get_triangular_iteration(45) == 10
    @test OscoNet.get_triangular_iteration(780.0) == 40
end
