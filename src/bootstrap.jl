"""
# Arguments

# Returns
"""
function bootstrap_hypothesis_test(
    data,
    n_bootstrap::Int;
    alpha=0.001
)
    n_genes, _ = size(data)

    _, cost_ng = find_best_psi_for_each_gene_pair(data)
    cost_permuted = get_permuted_cost(data, n_bootstrap)
    
    pvalues = get_pvalues(cost_ng,cost_permuted)
    pvalue_flatten = vec_triu_loop(pvalues)
    @assert length(pvalue_flatten) == n_genes * (n_genes-1) / 2 "Length of vector should be n_genes * (n_genes-1) / 2"
    
    qvalues_flatten, _ = qvalue_estimate(pvalue_flatten; verbose=true)   
    qvalues = get_symmetric_matrix_from_upper_triangular(size(pvalues)[1], qvalues_flatten)

    adjacency_matrix = qvalues .< alpha

    return adjacency_matrix, qvalues, cost_ng
end

"""

# Arguments

# Returns
"""
function find_best_psi_for_each_gene_pair(data)
    n_genes = size(data)[1]
    cost_ng = Array{Float64}(undef,n_genes,n_genes)
    psi_ng = Array{Float64}(undef,n_genes,n_genes)
    @inbounds @views begin
        for ix in 1:n_genes
            @floop for iy in (ix+1):n_genes
                res = find_minimum_distance(data[ix,:],data[iy,:])
                cost_ng[ix, iy] = res[1]
                psi_ng[ix, iy] = res[2]
            end
        end
    end
    return psi_ng, cost_ng
end


"""

# Arguments

- `data::`: G X N tensor of gene expression

- `n_permutations::Int`: number of bootstrap permutations

# Returns

- `cost_permuted::Array{Float64}`: permuted cost
"""
function get_permuted_cost(
    data,
    n_permutations::Int
)
    @assert n_permutations > 1 "Number of bootstrap permutations must be greater than 1"
    cost_permuted = Array{Float64}(
        undef,
        first(size(data)),
        first(size(data)),
        n_permutations
    )

    @inbounds @views begin
        for ix in 1:first(size(data))
            @floop for iy in (ix+1):first(size(data))
                for permutation_index in 1:n_permutations
                    cost_permuted[ix, iy, permutation_index] = first(find_minimum_distance(data[ix,:], shuffle(data[iy, :])))
                end
            end
        end
    end
    
    return cost_permuted
end

"""
A function to return the elements of the (strict) upper triangle of a matrix, M, as a vector
"""
function vec_triu_loop(M::AbstractMatrix{T}) where T
    m, n = size(M)
    @assert m == n
    l = n*(n+1) ÷ 2 - n
    v = Vector{T}(undef,l)
    k = 0
    @inbounds for i in 2:n
        for j in 1:i
            v[k + j] = M[j, i]
        end
        k += (i-1)
    end
    v
end

function get_symmetric_matrix_from_upper_triangular(n_genes::Int, flatten_vector)
    @assert length(flatten_vector) == n_genes*(n_genes-1)/2
    
    a = zeros(n_genes, n_genes)
    k = 1
    for j in 2:n_genes
        for i in 1:j-1
            a[i,j] = flatten_vector[k]
            k += 1
        end
    end
    
    a += a'
    a[diagind(a)] .= Inf  # put infinities on diagonal to protect caller from using these zero entries
    return a
end

