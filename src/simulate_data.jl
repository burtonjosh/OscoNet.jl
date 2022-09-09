"""
Generate synthetic data. Reproduced from the original [Oscope paper](https://www.nature.com/articles/nmeth.3549)

Included below is the original description from the supplementary material of that paper:

## Sim I: Oscope paper supplementary
1000 genes and 100 cells.  90 out of the 1000 genes were  simulated  as  oscillators.
The 90 oscillators were simulated in 3 frequency groups, each group contains 30 genes.
Group 1 and 3 follow the same order, while genes in group 2 follow another order.

The relative speeds of the 3 groups are proportional to 2:3:6.
Within each frequency group, genes were further simulated with strong and weak signals.
Half of the oscillatory genes were simulated as strong oscillators with `σ_g` = `σ_str`.
The other half were simulated as weak oscillators with `σ_g` = `2σ_str`.
Starting Ψ `Ψ_g` varies in different genes within a frequency group.
The remaining genes except the oscillators are called noise genes.
Noise genes were simulated as random Gaussian noise. The noise level was adjusted to
be comparable to the average noise signal among all oscillators.
The `σ_str` varies from 0.05 to 0.4 in 5 steps.

# Arguments

- `half_group::Int=15`: Half the size of each gene group. For example if `half_group=3`, each group will have 6 co-oscillating genes, 3 of which 
will be strong oscillators and 3 weak (double the noise)

- `n_genes::Int=1000`: Total number of genes

- `n_cells::Int=100`: Total number of cells

- `noise_level::Int=1`: Noise level index (1 to 6)

- `n_groups::Int=3`: Number of groups

# Returns

- `data::Matrix{<:AbstractFloat}`: Simulated gene expression data with shape `n_genes` x `n_cells` 

- `Ψ_g::Vector{<:AbstractFloat`: Vector of phases for each gene

- `ω::Vector{<:AbstractFloat}`: Vector of angular speeds for each gene
"""
function simulate_data(;
    half_group::Int = 15,
    n_genes::Int = 1000,
    n_cells::Int = 100,
    noise_level::Int = 1,
    n_groups::Int = 3,
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
    cell_names = ["C$i" for i = 0:n_cells-1]
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

    @assert !any(isnan, data)
    @assert length(gene_names) == n_genes
    @assert length(cell_names) == n_cells

    return data, Ψ_g, ω  # return GXN matrix

    # df = DataFrame(data, cell_names)
    # df[!,"gene_names"] = gene_names
    # return groupby(df,:gene_names), Ψ_g, ω
    # return df, Ψ_g, ω  # return GXN matrix

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