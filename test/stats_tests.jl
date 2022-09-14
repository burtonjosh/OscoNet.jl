using OscoNet, Random, LinearAlgebra

@testset "p-values and q-values" begin
    Random.seed!(1234)
    n_permutations = 100
    data, _, ω = simulate_data()
    _, cost_unpermuted = OscoNet.find_best_psi_for_each_gene_pair(data)
    cost_permuted = OscoNet.get_permuted_cost(data, n_permutations)

    pvalues = OscoNet.get_pvalues(cost_unpermuted, cost_permuted)

    @test size(pvalues) == (first(size(data)),first(size(data)))
    @test isapprox(sum(OscoNet.vec_triu_loop(pvalues) .== 0), 45; atol=2)

    pvalues_flatten = OscoNet.vec_triu_loop(pvalues)
    qvalues_flatten, π₀ = OscoNet.qvalue_estimate(pvalues_flatten)
    qvalues = OscoNet.construct_symmetric_matrix_from_upper_triangular_vector(qvalues_flatten)

    @test size(qvalues) == (first(size(data)),first(size(data)))
    @test issymmetric(qvalues)

    small_qvalues_flatten, small_π₀ = OscoNet.qvalue_estimate(pvalues_flatten[1:45])
    @test small_π₀ == 1

    small_qvalues_flatten, small_π₀ = OscoNet.qvalue_estimate(pvalues_flatten[1:45]; π₀=0.5)
    @test small_π₀ == 0.5

    adj_matrix_true = true_adj_matrix(ω)
    α_values = [1.0, 0.05, 0.01, 0.001]

    tpr, fpr, fdr = OscoNet.get_metrics_for_different_qvalue_thresholds(qvalues, adj_matrix_true, α_values)
    @test length(tpr) == length(α_values)
    @test length(fpr) == length(α_values)
    @test length(fdr) == length(α_values)

    @test first(tpr) == 1.0
    @test first(fpr) == 1.0
    @test first(fdr) == (190-45)/190

    # check fdr and fpr are zero if no gene pairs are discovered
    tpr, fpr, fdr = OscoNet.get_metrics_for_different_qvalue_thresholds(ones(20,20), adj_matrix_true, [0.01])
    @test first(fpr) == 0.0
    @test first(fdr) == 0.0
end