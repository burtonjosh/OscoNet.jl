using OscoNet

@testset "Distance function" begin
    xs = LinRange(-1, 1, 200)
    Ψ = π / 2
    X = sin.(xs)
    Y = sin.(xs .+ Ψ)

    @test isapprox(cooscillation_distance(X, Y, cos(Ψ)), 0.0; atol = 1e-15)
    @test cooscillation_distance(X, Y, cos(π)) ≈ 319.14974846118
end

@testset "Optimisation" begin
    xs = LinRange(-1, 1, 200)
    Ψ = π / 2
    X = sin.(xs)
    Y = sin.(xs .+ Ψ)

    minimum_distance, minimiser = find_minimum_distance(X, Y)

    @test isapprox(minimum_distance, 0.0; atol = 1e-15)
    @test minimiser ≈ π / 2
end
