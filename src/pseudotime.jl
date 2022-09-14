"""
# Arguments

- `data`

- `n_neighbors`

# Returns

- `pseudotime`

- `data_tsne_scaled`
"""
function estimate_pseudotime_using_tsne(data)
    data_tsne = predict(fit(TSNE, data))
    data_tsne_scaled = scale_data(data_tsne)

    # Fit circle
    centre = first(get_circle(data_tsne_scaled))
    data_tsne_scaled = data_tsne_scaled .- centre
    pseudotime_radians = atan.(data_tsne_scaled[1, :], data_tsne_scaled[2, :])
    
    # convert to [0, 1]
    pseudotime = (pseudotime_radians .- minimum(pseudotime_radians)) ./ (maximum(pseudotime_radians) - minimum(pseudotime_radians))
    return pseudotime, data_tsne_scaled
end

"""
Z-Score normalisation of a matrix of size (n_genes x n_cells), applied row-wise
"""
function scale_data(Y)
    μ = mean(Y, dims = 2)
    σ = std(Y, dims = 2)
    return (Y .- μ) ./ σ
end

"""
Calculate the distance of a set of 2D points from the centre (xc, yc) of a circle
"""
function distance_from_centre(data, centre)
    return sqrt.((data[1, :] .- centre[1]) .^ 2 .+ (data[2, :] .- centre[2]) .^ 2)
end

"""
Helper function for optimisation. Calculate the sum of the squared distances
between a set of data points and the circle centred at c=(xc, yc)
"""
function f_2(data, centre)
    radius = distance_from_centre(data, centre)
    return sum((radius .- mean(radius)) .^ 2)
end

"""
Calculate mean circle
"""
function get_circle(data)
    @assert size(data)[1] == 2 "Number of rows must be 2"

    res = optimize(x -> f_2(data, x), mean(data, dims = 2), BFGS(); autodiff = :forward)
    centre = Optim.minimizer(res)

    radius = distance_from_centre(data, centre)
    mean_radius = mean(radius)
    residual = sum((radius .- mean_radius) .^ 2)
    return centre, mean_radius, residual
end

"""
Evaluate roughness of pseudotime. This metric measures the smoothness of the gene expression profile
by looking at the differences between consecutive measurements. Smaller values indicate a smoother response. 
"""
function calc_roughness(x, pseudotime)
    i = sortperm(pseudotime)
    x = x[:, i]
    N = size(x)[2]
    S = std(x, dims = 2)
    return sqrt.((1 .+ sum((x[:, 1:(N-1)] - x[:, 2:N]) .^ 2, dims = 2)) ./ (N - 1)) ./ S
end

"""
Calculate time of peak expression for a specific gene given pseudotime
"""
function get_peak_time(data, clone_id, pseudotime; base_cycle_time = nothing)
    d = vec(Array(data[in(clone_id).(data.X0), :][!, 2:end]))
    if base_cycle_time == nothing
        base_cycle = d
    else
        base_cycle = d[base_cycle_time[1]:base_cycle_time[2]]
    end
    return pseudotime[argmax(base_cycle)]
end

"""
Calculate various metrics
"""
function calculate_metrics(data, comm_data, clone_id_list, pseudotime, realtime)
    peak_times = zeros(length(clone_id_list))
    roughness = zeros(length(clone_id_list))

    peak_times_true = zeros(length(clone_id_list))
    roughness_true = zeros(length(clone_id_list))

    for (id_index, id) in enumerate(clone_id_list)
        expression = comm_data[comm_data.clone_id.==id, :]
        peak_times[id_index] = get_peak_time(data, expression.clone_id, pseudotime)
        peak_times_true[id_index] = get_peak_time(data, expression.clone_id, realtime)
        roughness[id_index] =
            calc_roughness(vec(Array(data[in(expression.clone_id).(data.X0), :][:, 2:end]))', pseudotime)[1]
        roughness_true[id_index] = calc_roughness(
            vec(Array(data[in(expression.clone_id).(data.X0), :][:, 2:end]))',
            realtime,
        )[1]
    end

    # KS to test uniformity
    pseudotime_KS = HypothesisTests.ksstats(pseudotime, Uniform())[2]

    pseudotime_corr = corspearman(pseudotime, realtime)
    peak_corr = corspearman(peak_times, peak_times_true)

    return peak_times, peak_times_true, roughness, pseudotime_KS, pseudotime_corr, peak_corr
end
