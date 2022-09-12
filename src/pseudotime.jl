"""
# Arguments

- `data`

- `n_neighbors`

# Returns

- `pt`

- `X_init`
"""
function estimate_pseudotime_using_tsne(
    data
)
    X_init_embed = predict(fit(TSNE, data))
    X_init = scale_data(X_init_embed)  # zero mean, unit variance

    
    # Fit circle
    center, _, _ = get_circle(X_init)
    X_init = X_init .- center  # remove mean
    pt_radians = atan.(X_init[1, :], X_init[2, :])

    pt = (pt_radians .- minimum(pt_radians)) ./ (maximum(pt_radians) - minimum(pt_radians))  # convert to [0, 1]
    return pt, X_init
end

function scale_data(Y)
    μ = mean(Y,dims=2)
    σ = std(Y,dims=2)
    return (Y .- μ) ./ σ
end

""" calculate the distance of each 2D points from the center (xc, yc) """
function calc_R(xc, yc, dataC)
    return sqrt.((dataC[1, :] .- xc).^2 .+ (dataC[2, :] .- yc).^2)
end

""" calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
function f_2(c, dataC)
    Ri = calc_R(c[1], c[2], dataC)
    return sum((Ri .- mean(Ri)).^2)
end

"""
Calculate circle using lease circles
"""
function get_circle(dataC)
    @assert size(dataC)[1] == 2 "Number of rows must be 2"
    
    res = optimize(x->f_2(x,dataC), mean(dataC,dims=2), BFGS(); autodiff = :forward) # TODO
    # res = optimize(x->f_2(x,dataC), mean(dataC,dims=2))
    center_2 = Optim.minimizer(res)

    Ri_2 = calc_R(center_2..., dataC)
    R_2 = mean(Ri_2)
    residu_2 = sum((Ri_2 .- R_2).^2)
    return center_2, R_2, residu_2
end

"""
peak time function
"""
function get_peak_time(data, cloneid, pt; base_cycle_time=nothing)
    d = vec(Array(data[in(cloneid).(data.X0),:][!,2:end]))
    if base_cycle_time == nothing
        base_cycle = d
    else
        base_cycle = d[base_cycle_time[1]:base_cycle_time[2]]
    end
    return pt[argmax(base_cycle)]
end

"""
evaluate roughness of pseudotime
This metric measures the smoothness of the gene expression profile by looking at the differences
of consecutive measurements.
Smaller values indicate a smoother response. 
"""
function calc_roughness(x, pt)
    i=sortperm(pt)
    x = x[:, i]
    N = size(x)[2]
    S = std(x,dims=2)
    return sqrt.((1 .+ sum((x[:,1:(N-1)] - x[:,2:N]).^2, dims=2) ) ./ (N-1)) ./ S
end

"""
metrics function
"""
function calculate_metrics(data, commData, cloneidlistEval, pt, ptTruePeriod)
    # Measure fit
    peakTimes = zeros(length(cloneidlistEval))
    roughness = zeros(length(cloneidlistEval))

    peakTimes_true = zeros(length(cloneidlistEval))
    roughness_true = zeros(length(cloneidlistEval))

    for (ic, c) in enumerate(cloneidlistEval)
        a = commData[commData.CloneID .== c,:]
        peakTimes[ic] = get_peak_time(data, a.CloneID, pt)
        peakTimes_true[ic] = get_peak_time(data, a.CloneID, ptTruePeriod)
        roughness[ic] = calc_roughness(
            vec(Array(data[in(a.CloneID).(data.X0),:][:,2:end]))',
            pt
        )[1]
        roughness_true[ic] = calc_roughness(
            vec(Array(data[in(a.CloneID).(data.X0),:][:,2:end]))',
            ptTruePeriod
        )[1]
    end
    
    # KS test to uniformity
    ptKS = HypothesisTests.ksstats(pt,Uniform())[2]

    pseudoCorr = corspearman(pt, ptTruePeriod)
    peakCorr = corspearman(peakTimes, peakTimes_true)

    println("peak R=$peakCorr")
    println("pseudotime: true roughness = $(median(roughness .- median(roughness_true)))")
    println("corr pseudotime = $(pseudoCorr), ks=$ptKS")

    return peakTimes, peakTimes_true, roughness, ptKS, pseudoCorr, peakCorr
end