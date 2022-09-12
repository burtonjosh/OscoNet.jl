"""
Get pvalues from bootstrap

# Arguments

- `cost_unpermuted`: G X G - only upper triangular used
- `cost_permuted`: G X G X n permutations

#Returns

- `pvalues`
"""
function get_pvalues(cost_unpermuted, cost_permuted)

    n_genes = size(cost_unpermuted)[1]
    pvalues = fill(Inf, n_genes, n_genes)

    @views @inbounds begin
        @floop for ix in 1:n_genes
            for iy in (ix+1):n_genes
                pvalues[ix, iy] = percentilerank(cost_permuted[ix,iy, :], cost_unpermuted[ix,iy])/100.0
            end
        end
    end

    return pvalues
end

"""
Estimates qvalues from pvalues

# Arguments

- `pvalues`: 

- `m`: number of tests. If not specified m = length(pvalues)

- `verbose::Bool`: print verbose messages (default is false)

- `pi0`: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
"""
function qvalue_estimate(pvalues; m=nothing, verbose=false, pi0=nothing)

    @assert minimum(pvalues) >= 0 && maximum(pvalues) <= 1 "p-values should be between 0 and 1"

    if m == nothing
        m = length(pvalues)
    end

    # if the number of hypotheses is small, just set pi0 to 1
    if length(pvalues) < 100 && pi0 == nothing
        pi0 = 1.0
    elseif pi0 != nothing
        pi0 = pi0
    else
        # evaluate pi0 for different lambdas
        lambdas = LinRange(0., 0.9, 91)
        # counts = [sum(pvalues .> i) for i in lambdas]
        pi0 = sum(pvalues .> lambdas[end]) / (m*(1-lambdas[end]))


        # not obvious to me why a cubic spline is needed here 
        # fit natural cubic spline
        # println(pi0[end])
        # pi0 = [counts[l] / (m*(1-lambdas[l])) for l in 1:length(lambdas)]
        # tck = cubic_spline_interpolation(lambdas, pi0)
        # println(tck(lambdas[end]))
        # # println(tck(1.0))
        # pi0 = tck(lambdas[end])
        # println(pi0)

        if verbose
            println("qvalues pi0=",round(pi0; digits=3),", estimated proportion of null features ")
        end

        if pi0 > 1
            pi0 = 1.0

            if verbose
                println("got pi0 > 1 (", round(pi0; digits=3),") while estimating qvalues, setting it to 1")
            end

        end
    end

    @assert pi0 >= 0 && pi0 <= 1 "pi0 is not between 0 and 1: $pi0"

    p_ordered = sortperm(pvalues)
    pvalues = pvalues[p_ordered]
    qv = pi0 * m/length(pvalues) .* pvalues
    qv[end] = min(qv[end], 1.0)

    for i in (length(pvalues)-1):-1:1
        qv[i] = min(pi0*m*pvalues[i]/(i), qv[i+1])
    end

    # reorder qvalues
    qv_temp = copy(qv)
    qv = zero(qv)
    qv[p_ordered] = qv_temp

    return qv, pi0
end

"""
Get true positive, false discovery and false positive rates
:param qvalues: qvalue G X G matrix where G=number of genes
:param adjMatrix_true: GXG binary adjacency matrix of true gene pairs
:param alpha_values: one dimensional vector for threshold values
:return: TPR, FDR, FPR vector of same size as alpha_values
"""
function get_metrics_for_different_qvalue_thresholds(
qvalues,
adjMatrix_true,
alpha_values
)

    # FDR and TPR for different thresholds
    G = size(adjMatrix_true)[1]

    TPR, FDR, FPR = zeros(length(alpha_values)), zeros(length(alpha_values)), zeros(length(alpha_values))
    
    for (iaT, aT) in enumerate(alpha_values)
        adjMatrixBootstrapQvalue = zeros(G, G)
        adjMatrixBootstrapQvalue[qvalues .< aT] .= 1

        # println(adjMatrixBootstrapQvalue)

        TPR[iaT], FPR[iaT], FDR[iaT], _, _, _, _ = calculate_error_rates(
            Bool.(adjMatrixBootstrapQvalue),
            adjMatrix_true
        )
    end
    return TPR, FPR, FDR
end

function calculate_error_rates(adjMatrixEstimated, adjMatrixTrue)
    
    TP = sum(adjMatrixEstimated .& adjMatrixTrue)  # true positive
    TN = sum(.!(adjMatrixEstimated) .& .!(adjMatrixTrue))
    FP = sum(adjMatrixEstimated .& .!(adjMatrixTrue))  # false positive
    FN = sum(.!(adjMatrixEstimated) .& adjMatrixTrue)
    # True positive rate -TP divided by number of true oscillating pairs
    TPR = TP / sum(adjMatrixTrue)
    # False discovery rate -FP divided by number of gene pairs reported
    FPR = FP / (FP+TN)
    
    if (sum(adjMatrixEstimated) == 0)
        # No co-osc genes found - so FDR is 0
        FDR = 0.
        FPR = 0.
    else
        FDR = FP / sum(adjMatrixEstimated)
    end
    
    @assert !any(isnan, TPR)#np.all(~np.isnan(TPR))
    @assert !any(isnan, FDR)#np.all(~np.isnan(FDR))
    
    return TPR, FPR, FDR, TP, TN, FP, FN

end