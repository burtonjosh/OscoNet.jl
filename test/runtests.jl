using SafeTestsets

@safetestset "Data simulation" begin
    include("data_simulation_tests.jl")
end

@safetestset "Oscope" begin
    include("oscope_tests.jl")
end

@safetestset "Bootstrap" begin
    include("bootstrap_tests.jl")
end

@safetestset "Edge network" begin
    include("network_tests.jl")
end

@safetestset "Statistics" begin
    include("stats_tests.jl")
end

@safetestset "Pseudotime" begin
    include("pseudotime_tests.jl")
end
