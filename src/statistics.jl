"""
Get p-values from the bootstrap hypothesis test

# Arguments

- `cost_unpermuted`: A matrix of size (n_genes X n_genes) whose upper triangular `(x,y)`-th element defines the
distance between gene `x` and gene `y`

- `cost_permuted`: An array of size (n_genes X n_genes X n_permutations) whose upper triangular `(x,y,n)`-th element
defines the distance between gene `x` and the `n`-th permutation of gene `y`

# Returns

- `pvalues`: A matrix of size (n_genes X n_genes)

"""
function get_pvalues(cost_unpermuted, cost_permuted)

    n_genes = size(cost_unpermuted)[1]
    pvalues = fill(Inf, n_genes, n_genes)

    @views @inbounds begin
        @floop for ix = 1:n_genes
            for iy = (ix+1):n_genes
                pvalues[ix, iy] =
                    percentilerank(cost_permuted[ix, iy, :], cost_unpermuted[ix, iy]) /
                    100.0
            end
        end
    end

    return pvalues
end

"""
Estimates q-values from p-values, see [https://github.com/nfusi/qvalue](https://github.com/nfusi/qvalue)

# Arguments

- `pvalues`: A vector of p-values

- `n_tests`: number of tests. If not specified n_tests = length(pvalues)

- `π₀`: Estimate of the proportion of features that are truly null. If not specified,
it's estimated as suggested in Storey and Tibshirani, 2003.

# Returns

- `qvalues`

- `π₀`
"""
function qvalue_estimate(pvalues; n_tests = nothing, π₀ = nothing)

    @assert minimum(pvalues) >= 0 && maximum(pvalues) <= 1 "p-values should be between 0 and 1"

    if n_tests == nothing
        n_tests = length(pvalues)
    end

    # if the number of hypotheses is small, just set π₀ to 1
    if length(pvalues) < 100 && π₀ == nothing
        π₀ = 1.0
    elseif π₀ != nothing
        π₀ = π₀
    else
        # evaluate π₀ for different lambdas
        lambdas = LinRange(0.0, 0.9, 91)
        # counts = [sum(pvalues .> i) for i in lambdas]
        π₀ = sum(pvalues .> lambdas[end]) / (n_tests * (1 - lambdas[end]))


        # not obvious to me why a cubic spline is needed here 
        # fit natural cubic spline
        # println(π₀[end])
        # π₀ = [counts[l] / (m*(1-lambdas[l])) for l in 1:length(lambdas)]
        # tck = cubic_spline_interpolation(lambdas, π₀)
        # println(tck(lambdas[end]))
        # # println(tck(1.0))
        # π₀ = tck(lambdas[end])
        # println(π₀)

        if π₀ > 1
            π₀ = 1.0
        end
    end

    @assert π₀ >= 0 && π₀ <= 1 "π₀ is not between 0 and 1: $π₀"

    p_ordered = sortperm(pvalues)
    pvalues = pvalues[p_ordered]
    qvalues = π₀ * n_tests / length(pvalues) .* pvalues
    qvalues[end] = min(qvalues[end], 1.0)

    for i = (length(pvalues)-1):-1:1
        qvalues[i] = min(π₀ * n_tests * pvalues[i] / (i), qvalues[i+1])
    end

    # reorder qvalues
    qv_temp = copy(qvalues)
    qvalues = zero(qvalues)
    qvalues[p_ordered] = qv_temp

    return qvalues, π₀
end

"""
Get true positive, false positive and false discovery rates

# Arguments

- `qvalues`: q-value matrix (n_genes X n_genes)

- `adj_matrix_true`: binary adjacency matrix of true gene pairs (n_genes X n_genes)

- `α_values`: one dimensional vector for threshold values

# Returns

- `true_positive_rate`:

- `false_positive_rate`:

- `false_discovery_rate`:
"""
function get_metrics_for_different_qvalue_thresholds(qvalues, adj_matrix_true, α_values)

    # false_discovery_rate and true_positive_rate for different thresholds
    n_genes = first(size(adj_matrix_true))

    true_positive_rate = zeros(length(α_values))
    false_positive_rate = zeros(length(α_values))
    false_discovery_rate = zeros(length(α_values))

    for (α_index, α) in enumerate(α_values)
        adj_matrix_qvalue = zeros(n_genes, n_genes)
        adj_matrix_qvalue[qvalues.<α] .= 1

        # calculate error rates on upper triangle only using `vec_triu_loop` to avoid diagonal and duplicating pairs
        true_positive_rate[α_index],
        false_positive_rate[α_index],
        false_discovery_rate[α_index],
        _,
        _,
        _,
        _ = calculate_error_rates(vec_triu_loop(Bool.(adj_matrix_qvalue)), vec_triu_loop(adj_matrix_true))
    end
    return true_positive_rate, false_positive_rate, false_discovery_rate
end

"""
# Arguments

-

# Returns

-
"""
function calculate_error_rates(adj_matrix_estimate, adj_matrix_true)

    true_positive = sum(adj_matrix_estimate .& adj_matrix_true)
    true_negative = sum(.!(adj_matrix_estimate) .& .!(adj_matrix_true))
    false_positive = sum(adj_matrix_estimate .& .!(adj_matrix_true))
    false_negative = sum(.!(adj_matrix_estimate) .& adj_matrix_true)

    true_positive_rate = true_positive / sum(adj_matrix_true)
    false_positive_rate = false_positive / (false_positive + true_negative)

    if sum(adj_matrix_estimate) == 0
        false_discovery_rate = 0.0
        false_positive_rate = 0.0
    else
        false_discovery_rate = false_positive / sum(adj_matrix_estimate)
    end

    @assert !any(isnan, true_positive_rate)
    @assert !any(isnan, false_discovery_rate)

    return true_positive_rate,
    false_positive_rate,
    false_discovery_rate,
    true_positive,
    true_negative,
    false_positive,
    false_negative

end
