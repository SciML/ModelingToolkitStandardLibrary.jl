"""
2D rigid body with rectangular shape. The parameters `s_start` and `s_end` are two points in 2D space. Together with the parameter `width`, the vector drawn between these two points defines the spatial extents of our rigid body. This "two-point" vector implementation differs from Modelica's implementation of the same component. 
In the Modelica implementation, one point is used to define the position of the origin of a new frame of reference, and a second point is used to create a vector from the origin of this new frame.

# Parameters:
- `s_start`: [m] Point in 2D space that defines the starting position of our rigid body
- `s_end`: [m] Point in 2D space that defines the ending position of our rigid body
- `w_start`: [rad/s] Initial angular velocity of the ends of our rigid body, to be defined in a direction perpendicular to the vector that forms the body itself
- `a_start`: [rad/s^2] Initial angular acceleration of the ends of our rigid body, to be defined in a direction perpendicular to the vector that forms the body itself
- `width`: [m] Width of our rigid body
- `ρ`: [kg/m^3] Material density of our rigid body

"""
struct Box
    s_start::Vector{Float64} = [0.0, 0.0]
    s_end::Vector{Float64} = [1.0, 1.0]
    w_start::Float64 = 0.0
    a_start::Float64 = 0.0
    width::Float64 = 1.0
    ρ::Float64 = 2710.0

    @parameters t ρ = ρ
    @parameters length = sqrt ( ( s_start[1] - s_end[1] ) ^ 2 + ( s_start[2] - s_end[2] ) ^ 2 )
    @parameters mass = ρ * length * width

    D = Differential(t)

    sts = @variables begin
        x(t) = ( s_start[1] + s_end[1] ) / 2
        y(t) = ( s_start[2] + s_end[2] ) / 2
        w(t) = w_start
        a(t) = a_start
    end

    eqs = [
        D(s) ~ w
        D(w) ~ a
        m * a ~ 0
    ]
end

# ODESystem(eqs, t, sts, [m]; name=name)