using LinearAlgebra

"""
Uniform gravitational field, defined here as a constant field vector `g`. This vector is produced from two inputs: field magnitude and field 
direction (defaulted to 9.81 and [0, -1], respectively. All units are in S.I. units). Furthermore, a single invariant is enforced with an inner 
constructor method — that is, the provided field magnitude must equal or exceed 0.0. Note that the provided direction vector is normalized 
within our constructor, and so it is not necessary to input a unit vector.

# Parameters:
- `g`: [m/s²] Constant gravitational field vector
"""
struct UniformField
    g::Vector{Float64}
    UniformField(g::Float64 = 9.81, n::Vector{Float64} = [0.0,-1.0]) = g >= 0.0 ? new(g * normalize(n)) : error("Magnitude of gravitational acceleration must exceed zero")
end