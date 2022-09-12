"""
Generate synthetic data. Reproduced from the original [Oscope paper](https://www.nature.com/articles/nmeth.3549)

This function can be used to simulate fake gene expression data with which to test the other functions.

# Arguments

- `half_group::Int=10`: Half the size of each gene group. For example if `half_group=3`, each group will have 6 co-oscillating genes, 3 of which 
will be strong oscillators and 3 weak (double the noise)

- `n_genes::Int=80`: Total number of genes

- `n_cells::Int=1000`: Total number of cells

- `noise_level::Int=1`: Noise level index (1 to 6)

- `n_groups::Int=2`: Number of groups

# Returns

- `data::Matrix{<:AbstractFloat}`: Simulated gene expression data with shape `n_genes` x `n_cells` 

- `Ψ_g::Vector{<:AbstractFloat`: Vector of phases for each gene

- `ω::Vector{<:AbstractFloat}`: Vector of angular speeds for each gene
"""
function simulate_data(;
    half_group::Int = 10,
    n_genes::Int = 80,
    n_cells::Int = 1000,
    noise_level::Int = 1,
    n_groups::Int = 2,
)
    @assert n_groups <= 3

    # Construct oscillatory groups
    σ_str_level = [0.05, 0.1, 0.2, 0.3, 0.4, 0.6]
    σ_str = σ_str_level[noise_level]

    # two different orders
    t1 = LinRange(0, 2pi, n_cells)
    t2 = shuffle(t1)

    data = fill(NaN, n_genes, n_cells)

    # genes per weak/strong oscillatory group
    # Group 1
    # cell_names = ["C$i" for i = 0:n_cells-1] # TODO make use of cell names - dataframe?
    gene_names = vcat(["G1S0$i" for i = 0:half_group-1], ["G1W0$i" for i = half_group:2half_group-1])

    Ψ_g = zeros(n_genes)
    ω = zeros(n_genes)

    white_noise_index = 2half_group + 1

    for i = 1:half_group  # strong oscillators
        starting_Ψ = 2pi * rand()
        Ψ_g[i] = starting_Ψ
        ω[i] = 2
        data[i, :] = sin.(2t1 .+ starting_Ψ) .+ σ_str * randn(n_cells)
    end
    for i = half_group+1:2half_group  # weak oscillators
        starting_Ψ = 2pi * rand()
        Ψ_g[i] = starting_Ψ
        ω[i] = 2
        data[i, :] = sin.(2t1 .+ starting_Ψ) .+ 2σ_str * randn(n_cells)
    end

    if n_groups >= 2
        gene_names =
            vcat(gene_names, ["G2S0$i" for i = 2half_group:3half_group-1], ["G2W0$i" for i = 3half_group:4half_group-1])
        white_noise_index = 4half_group + 1

        for i = (2half_group+1):3half_group  # strong oscillators
            starting_Ψ = 2pi * rand()
            Ψ_g[i] = starting_Ψ
            ω[i] = 3
            data[i, :] = sin.(3t2 .+ starting_Ψ) .+ σ_str * randn(n_cells)
        end
        for i = (3half_group+1):4half_group  # weak oscillators
            starting_Ψ = 2pi * rand()
            Ψ_g[i] = starting_Ψ
            ω[i] = 3
            data[i, :] = sin.(3t2 .+ starting_Ψ) .+ 2σ_str * randn(n_cells)
        end
    end

    if n_groups >= 3
        # Group 3
        gene_names =
            vcat(gene_names, ["G3S0$i" for i = 3half_group:4half_group-1], ["G3W0$i" for i = 4half_group:5half_group-1])
        white_noise_index = 6half_group + 1

        for i in (4half_group+1:5half_group)  # strong oscillators
            starting_Ψ = 2pi * rand()
            Ψ_g[i] = starting_Ψ
            ω[i] = 6
            data[i, :] = sin.(6t1 .+ starting_Ψ) .+ σ_str * randn(n_cells)
        end
        for i = (5half_group+1):6half_group  # weak oscillators
            starting_Ψ = 2pi * rand()
            Ψ_g[i] = starting_Ψ
            ω[i] = 6
            data[i, :] = sin.(6t1 .+ starting_Ψ) .+ 2σ_str * randn(n_cells)
        end
    end
    
    # white noise genes
    gene_names = vcat(gene_names, ["R$i" for i = white_noise_index:n_genes])
    for w = white_noise_index:n_genes  # use i index from above where it stopped
        Ψ_g[w] = NaN
        ω[w] = NaN
        data[w, :] = maximum([3 / 2 * σ_str, 1]) * randn(n_cells)
    end

    return data, Ψ_g, ω
end


"""
Construct the ground truth adjacency matrix for a given numer of genes and
`ω` vector. Intended for use with output from the [`simulate_data`](@ref simulate_data) function.

# Arguments

- `n_genes::Integer`: Total number of genes

- `ω::Vector{<:AbstractFloat}`: Vector of angular speeds

# Returns

- `adj_mat::Matrix{Bool}`: Boolean matrix which is 1 if genes co-oscillate, and 0 otherwise
"""
function true_adj_matrix(n_genes::Int, ω::Vector{<:AbstractFloat})
    adj_mat = fill(false, n_genes, n_genes)
    @views @inbounds begin
        for i = 1:n_genes
            for j = i+1:n_genes
                adj_mat[i, j] = (ω[i] == ω[j])
                adj_mat[j, i] = (ω[i] == ω[j])
            end
        end
    end
    return adj_mat
end