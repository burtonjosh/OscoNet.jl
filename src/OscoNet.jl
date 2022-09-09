module OscoNet

using DataFrames, FLoops, LinearAlgebra, Optim, Random, StatsBase

include("oscope.jl")
include("bootstrap.jl")
include("simulate_data.jl")
include("edge_network.jl")
include("statistics.jl")

export simulate_data, true_adj_matrix
export cooscillation_distance, find_minimum_distance
export bootstrap_hypothesis_test
export create_edge_network

end
