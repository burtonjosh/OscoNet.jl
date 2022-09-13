"""
Perform a non-parameteric bootstrap hypothesis test on gene expression data to 
identify co-oscillating gene pairs.

# Arguments

- `data::Matrix`: Matrix of size (n_genes x n_cells), whose elements represent gene expression.

- `n_permutations::Int`: number of bootstrap permutations

- `α::AbstractFloat`

# Returns

- `adjacency_matrix::`:

- `qvalues::`:

- `cost::`:
"""
function bootstrap_hypothesis_test(data, n_permutations::Int; α = 0.001::AbstractFloat)
    _, cost = find_best_psi_for_each_gene_pair(data)
    cost_permuted = get_permuted_cost(data, n_permutations)

    pvalues = get_pvalues(cost, cost_permuted)
    pvalue_flatten = vec_triu_loop(pvalues)

    qvalues_flatten, _ = qvalue_estimate(pvalue_flatten)
    qvalues = construct_symmetric_matrix_from_upper_triangular_vector(qvalues_flatten)

    adjacency_matrix = qvalues .< α

    return adjacency_matrix, qvalues, cost
end

"""
Compute the minimum distance between gene `x` and gene `y`, for
a set of genes `x, y ∈ X`. 

# Arguments

- `data::Matrix`: Matrix of size (n_genes x n_cells), whose elements represent gene expression.

# Returns

- `Ψ::Array{Float64}`: Array of size (n_genes x n_genes), whose `(x,y)`-th element
represents the phase shift Ψ between gene `x` and gene `y`.

- `distance::Array{Float64}`: Array of size (n_genes x n_genes), whose `(x,y)`-th element
represents the minimum distance `d` between gene `x` and gene `y`.
"""
function find_best_psi_for_each_gene_pair(data)
    n_genes = size(data)[1]
    cost = Array{Float64}(undef, n_genes, n_genes)
    Ψ = Array{Float64}(undef, n_genes, n_genes)
    @inbounds @views begin
        for ix = 1:n_genes
            @floop for iy = (ix+1):n_genes
                res = find_minimum_distance(data[ix, :], data[iy, :])
                cost[ix, iy] = res[1]
                Ψ[ix, iy] = res[2]
            end
        end
    end
    return Ψ, cost
end


"""
Compute the minimum distance between gene `x` and the `n`-th random permutation of gene `y`, for
a set of genes `x, y ∈ X`. 

# Arguments

- `data::Matrix`: Matrix of size (n_genes x n_cells), whose elements represent gene expression.

- `n_permutations::Int`: number of bootstrap permutations

# Returns

- `cost_permuted::Array{Float64}`: Array of size (n_genes x n_genes x n_permutations), whose `(x,y,n)`-th element
represents the minimum distance between gene `x` and the `n`-th random permutation of gene `y`.
"""
function get_permuted_cost(data, n_permutations::Int)
    @assert n_permutations > 1 "Number of bootstrap permutations must be greater than 1"
    cost_permuted =
        Array{Float64}(undef, first(size(data)), first(size(data)), n_permutations)

    @inbounds @views begin
        for ix = 1:first(size(data))
            @floop for iy = (ix+1):first(size(data))
                for permutation_index = 1:n_permutations
                    cost_permuted[ix, iy, permutation_index] =
                        first(find_minimum_distance(data[ix, :], shuffle(data[iy, :])))
                end
            end
        end
    end

    return cost_permuted
end

"""
A function to return the elements of the (strict) upper triangle of a matrix, M, as a vector.

# Examples
```jldoctest
julia> using OscoNet

julia> A = [1 2 3; 4 5 6; 7 8 9]
3×3 Matrix{Int64}:
 1  2  3
 4  5  6
 7  8  9

julia> OscoNet.vec_triu_loop(A)
3-element Vector{Int64}:
 2
 3
 6
```
"""
function vec_triu_loop(M::AbstractMatrix{T}) where {T}
    m, n = size(M)
    @assert m == n
    l = n * (n + 1) ÷ 2 - n
    v = Vector{T}(undef, l)
    vector_index = 1
    @inbounds for i = 2:n
        for j = 1:(i-1)
            v[vector_index] = M[j, i]
            vector_index += 1
        end
    end
    v
end

"""
Construct a symmetric matrix with Inf diagonal values using a vector whose length is a triangular
number. This vector is intended to be the output of `vec_triu_loop`.

# Examples
```jldoctest
julia> using OscoNet

julia> A = [1 2 3; 4 5 6; 7 8 9]
3×3 Matrix{Int64}:
 1  2  3
 4  5  6
 7  8  9

julia> A_flatten = OscoNet.vec_triu_loop(A)
3-element Vector{Int64}:
 2
 3
 6

julia> OscoNet.construct_symmetric_matrix_from_upper_triangular_vector(A_flatten)
3×3 Matrix{Float64}:
 Inf    2.0   3.0
  2.0  Inf    6.0
  3.0   6.0  Inf
```
"""
function construct_symmetric_matrix_from_upper_triangular_vector(flatten_vector)
    @assert istriangular(length(flatten_vector)) "length of input vector is not a triangular number"

    n_rows = get_triangular_iteration(length(flatten_vector))
    a = zeros(n_rows, n_rows)
    k = 1
    for j = 2:n_rows
        for i = 1:j-1
            a[i, j] = flatten_vector[k]
            k += 1
        end
    end

    a += a'
    a[diagind(a)] .= Inf  # put infinities on diagonal to protect caller from using these zero entries

    return a
end

function get_triangular_iteration(x::Number)
    @assert istriangular(x) "number is not triangular"
    (Int(sqrt(8x + 1)) + 1) ÷ 2
end

"""
    istriangular(x::Number) -> Bool

Test whether `x` is a triangular number.

By solving the quadratic formula for the the `m`-th triangular number `n = m(m-1)/2`,
it can be shown that an integer is triangular if and only if `8n+1` is a perfect square.
In other words, if and only if `sqrt(8n+1)` is an integer.

# Examples
```jldoctest
julia> using OscoNet

julia> OscoNet.istriangular(45.0)
true
```
"""
function istriangular(x::Number)
    @assert x >= 0 "number must be >= 0"
    isinteger(sqrt(8x + 1))
end
