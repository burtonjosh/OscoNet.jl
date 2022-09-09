"""
Distance function for two co-oscillating genes across `N` cells. We assume genes have identical
sinusoidal expression profiles, except for a phase shift Ψ.

# Arguments

- `X`: Expression of gene X from `N` cells
- `Y`: Expression of gene Y from the same `N` cells
- `cosΨ::AbstractFloat`: Cosine of the phase shift, Ψ

# Returns

- `d::AbstractFloat`: Distance between the two genes for a given Ψ
"""
function cooscillation_distance(
    X,
    Y,
    cosΨ::AbstractFloat,
)
    sum(@. (X^2 + Y^2 - 2cosΨ * X * Y - 1 + cosΨ^2)^2)
end


"""
Find the minimum distance between two genes X and Y,
and the value of Ψ that minimises it.

# Arguments

- `X`: Expression of gene X from `N` cells
- `Y`: Expression of gene Y from the same `N` cells

# Returns

- `minimum::AbstractFloat`: Minimum distance between gene X and Y
- `minimiser::AbstractFloat`: Value of Ψ that minimises the distance
"""
function find_minimum_distance(X, Y)
    res = optimize(x->cooscillation_distance(X, Y, x), -1.0, 1.0)
    return Optim.minimum(res), acos(Optim.minimizer(res))
end