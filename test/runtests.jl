using SafeTestsets

@safetestset "Data simulation" begin
    include("data_simulation_tests.jl")
end
@safetestset "Data simulation" begin
    include("oscope_tests.jl")
end
@safetestset "Data simulation" begin
    include("bootstrap_tests.jl")
end
