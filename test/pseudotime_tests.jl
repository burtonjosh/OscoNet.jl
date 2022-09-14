using OscoNet, Random, StatsBase

@testset "Pseudotime" begin
    Random.seed!(1234)
    n_cells=100
    data, _, _ = simulate_data(;n_cells=n_cells)

    pseudotime, data_tsne_scaled = estimate_pseudotime_using_tsne(data)

    @test length(pseudotime) == n_cells
    @test size(data_tsne_scaled) == (2,n_cells)
    @test all(pseudotime .<= 1) && all(pseudotime .>= 0)

end

@testset "Data manipulation" begin
    Random.seed!(1234)
    data, _, _ = simulate_data()
    
    scaled_data = OscoNet.scale_data(data)
    @test all(isapprox.(mean(scaled_data,dims=2), 0; atol=1e-15))
    @test all(isapprox.(std(scaled_data,dims=2), 1))

    unit_circle_points = hcat(
        LinRange(0,1,10),
        [√(1-x^2) for x in LinRange(0,1,10)]
    )'
    centre= [0.,0.]

    @test OscoNet.distance_from_centre(unit_circle_points, centre) == ones(10)
    @test OscoNet.f_2(unit_circle_points, centre) == 0

    centre_guess, radius_guess, residual = OscoNet.get_circle(unit_circle_points)
    @test isapprox(centre_guess, zeros(2); atol=1e-7)
    @test radius_guess ≈ 1.0
    @test isapprox(residual, 0.0; atol=1e-7)
end

@testset "Metrics" begin
    Random.seed!(1234)
    n_cells=100
    data, _, _ = simulate_data(;n_cells=n_cells)

    x = hcat(
        LinRange(0,1,100),
        LinRange(0,1,100)
    )'
    pseudotime = LinRange(0,1,100)

    @test all(OscoNet.calc_roughness(x, pseudotime) .== [0.34469099377285567, 0.34469099377285567])

    
end