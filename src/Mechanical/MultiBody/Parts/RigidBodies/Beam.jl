using ModelingToolkit, Parameters, LinearAlgebra

@inline function translate_vector_to_origin(P::Vector{Float64}, Q::Vector{Float64})
    return [Q[1]-P[1], Q[2]-P[2]]
end

@inline function compute_angle_between_vectors(P::Vector{Float64}, Q::Vector{Float64})
    return acos( dot(P,Q) / (norm(P) * norm(Q)) )
end

"""
2D rigid body with rectangular shape. The parameters `P` and `Q` are points in 2D space. Together with the parameter `width`, the vector drawn between points `P` and `Q` defines the spatial extents of our rigid body. This "two-point" vector implementation differs from Modelica's implementation of the same component. In the Modelica 
implementation, one point is used to define the position of the origin of a new frame of reference, and a second point is used to create a vector from the origin of this new frame.

# Parameters:
- `fixed`: [true/false] Boolean used to specify whether beam is fixed at a pivot point. If true, this pivot point is set at point `P` by default
- `P`: [m] Point in 2D space defining one end of our rigid body
- `Q`: [m] Point in 2D space defining opposite end of our rigid body
- `𝜏`: [N.m] Initial value of torque on our rigid body
- `ω`: [rad/s] Initial angular velocity of our rigid body
- `α`: [rad/s²] Initial angular acceleration of our rigid body
- `ρ`: [kg/m³] Material density of our rigid body
- `width`: [m] Width of our rigid body
- `length`: [m] Length of our rigid body, derived from points `P` and `Q`
- `m`: [kg] Mass of our rigid body, derived from `ρ`, `length`, and `width`
- `COM`: [m] Point in 2D space placed at the center of mass of our rigid body, derived from points `P` and `Q`
- `r`: [m] Distance between our pivot point and center of mass, derived from `length`
- `φ`: [rad] Initial value of absolute rotation angle of rigid body, measured from x-axis
- `θ`: [rad] Angle between the vector formed by points `P` and `Q`, and the gravity vector
"""
@with_kw struct Beam
    fixed::Bool = true
    P::Vector{Float64} = [0.0, 0.0]
    Q::Vector{Float64} = [1.0, 1.0]
    𝜏::Float64 = 0.0
    ω::Float64 = 0.0
    α::Float64 = 0.0
    ρ::Float64 = 2710.0
    width::Float64 = 1.0
    length::Float64 = norm(P-Q)
    m::Float64 = ρ * length * width 
    COM::Vector{Float64} = [ ( P[1] + Q[1] ) / 2, ( P[2] + Q[2] ) / 2 ]
    r::Float64 = length / 2
    φ::Float64 = compute_angle_between_vectors(translate_vector_to_origin(P,Q), [1.0,0.0])
    θ::Float64 = compute_angle_between_vectors(translate_vector_to_origin(P,Q), [0.0,-1.0])
end

beam = Beam()

@parameters t m=beam.m g=9.81 r=beam.r

sts = @variables begin
    P(t) = beam.P
    Q(t) = beam.Q
    𝜏(t) = beam.𝜏
    ω(t) = beam.ω
    α(t) = beam.α
    COM(t) = beam.COM
    φ(t) = beam.φ
    θ(t) = beam.θ
end

δ = Differential(t)

eqs = [δ(φ) ~ ω
       δ(ω) ~ α
       𝜏 ~ - m * g * r * sin(θ)
       δ(P) ~ 0]

@named de = ODESystem(eqs, t, sts, [m, g, r])