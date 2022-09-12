"""
Construct an undirected graph between gene pairs, to be used for identifying
'clusters' or communities of co-oscillating genes.

# Arguments

- `bootstrap_adj_mat::`:

- `weight_matrix::`:

- `gene_names::`:

# Returns

- `edge_network::DataFrame`: 

"""
function create_edge_network(
    bootstrap_adj_mat,
    weight_matrix,
    gene_names
)
    @assert size(bootstrap_adj_mat) == size(weight_matrix)
    @assert isapprox(bootstrap_adj_mat, bootstrap_adj_mat') "not symmetric"

    n_genes = size(weight_matrix)[1]
    @assert length(gene_names) == n_genes

    # Create edge representation
    edge_network = DataFrame(gene_1=String[], gene_2=String[], weight=Float64[])

    for i in 1:(n_genes-1)
        for j in (i+1):n_genes
            if bootstrap_adj_mat[i, j] == 1
                push!(edge_network,[gene_names[i], gene_names[j], weight_matrix[i, j]])
            end
        end
    end
    return edge_network
end