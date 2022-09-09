"""
CreateEdgeNetwork - Create Edge file.
This is needed before hypothesis test q-value derived adjacency matrix
can be consumed by R network analysis code.
Return a pandas dataframe with 3 columns, two gene names for the gene-pair and the cost value
"""
function create_edge_network(
    bootstrap_adj_mat,
    weight_matrix,
    gene_names
)
    @assert size(bootstrap_adj_mat) == size(weight_matrix)
    # we remove significant pairs that are not symmetric
    @assert isapprox(bootstrap_adj_mat, bootstrap_adj_mat') "not symmetric"
    G = size(weight_matrix)[1]
    nt = G*(G-1)  # number of tests without diagonal
    # println("Sparseness $(sum(bootstrap_adj_mat)/nt)") 
    # Get gene names
    @assert length(gene_names) == G

    # Create edge representation
    nSignificantPairs = sum(bootstrap_adj_mat) / 2  # symmetric matrix
    @assert isinteger(nSignificantPairs)

    edgeNetwork = DataFrame(gene1=String[], gene2=String[], weight=Float64[]) # np.empty((int(nSignificantPairs), 3), dtype='string, string, float64')

    for i in 1:(G-1)
        for j in (i+1):G
            if bootstrap_adj_mat[i, j] == 1
                push!(edgeNetwork,[gene_names[i], gene_names[j], weight_matrix[i, j]])
            end
        end
    end
    return edgeNetwork
end