"""
2D rigid body with rectangular shape. The parameters `P` and `Q` are two points in 2D space. Together with the parameter `width`, the vector drawn between points `P` and `Q` defines the spatial extents of our rigid body. This "two-point" vector implementation differs from Modelica's implementation of the same component. In the Modelica 
implementation, one point is used to define the position of the origin of a new frame of reference, and a second point is used to create a vector from the origin of this new frame.

# Parameters:
- `P`: [m] Point in 2D space defining one end of our rigid body
- `Q`: [m] Point in 2D space defining opposite end of our rigid body
- `ω`: [rad/s] Initial angular velocity of our rigid body. `ω` is to be defined as a vector, placed at opposite ends of our rigid body. This forms a pair of vectors that are anti-parallel
- `α`: [rad/s²] Initial angular acceleration of our rigid body. `α` is to be defined as a vector, placed at opposite ends of our rigid body. This forms a pair of vectors that are anti-parallel
- `width`: [m] Width of our rigid body
- `ρ`: [kg/m³] Material density of our rigid body
"""
using LinearAlgebra, SymPy

struct Box
    P::Vector{Float64} = [0.0, 0.0]
    Q::Vector{Float64} = [1.0, 1.0]
    ω::Float64 = 0.0
    α::Float64 = 0.0
    width::Float64 = 1.0
    ρ::Float64 = 2710.0

    @parameters t ρ = ρ
    @parameters length = norm(P-Q)
    @parameters COM = [P[1] + Q[1] / 2, P[2] + Q[2] / 2]
    @parameters mass = ρ * length * width

    D = Differential(t)

    sts = @variables begin
        P(t) = P
        Q(t) = Q
        ω(t) = ω
        α(t) = α
    end
end

@inline function compute_unit_normal_vector(P::Vector{Float64}, Q::Vector{Float64})
    Q = Q - P
    y = symbols("y")
    R = [1.0, y]
    y = solveset(Eq(dot(Q,R), 0), y)
    R[2] = y.args[1]
    return normalize([R[1], R[2]])
end

@inline function translate_vector(vector::Vector{Float64}, point::Vector{Float64})
    P = [0.0, 0.0] + point
    Q = vector + point
    return P, Q
end