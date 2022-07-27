include("GravityFields.jl")

"""
2D Cartesian coordinate system. Multiple instances of `Cartesian2D` are related via a tree structure, hence the existence of the parameters `parent` 
and `children`. When modeling a mechanical system, it often becomes necessary to first define a world that acts as absolute reference for our measurements, 
component placements, etc. In this program, a world is created from the default values for all our parameters — that is, `world = Cartesian2D()`. A world 
frame possesses no parent, initially has no children, and cannot be transformed or initialized with properties other than a gravitational field (hence, 
`translation`, `rotation`, `v`, etc. parameters are all defaulted to `nothing`). Again, in this program, world frames are defaulted to possess a uniform
gravitational field, defined in `GravityFields.jl` as a constant field vector `g` with magnitude 9.81 m/s² and direction vector [0, -1] (i.e. the g-field
is defaulted to act along the negative y-axis). 

Child frames can be added to the world by passing `world` as parent, and specifying whatever transformations and properties that are needed — all 
transformations are taken as being referenced from the origin of the world frame. For instance, `frame_a = Cartesian2D(parent=world, transformation=[1.0, 1.0])` 
creates a reference frame whose origin is placed at `(x, y) = (1.0, 1.0)` of the world frame. Importantly, when instantiating `frame_a`, `g` is inherited from 
`world` and `frame_a` is pushed onto the`children` stack of `world`.

# Parameters:
- `parent`: Parent frame of reference frame
- `children`: Child frame(s) of reference frame
- `translation`: [m] Translation vector of reference frame origin with respect to the parent
- `rotation`: [rad] Rotation of reference frame with respect to the parent
- `v`: [m/s] Velocity of reference frame with respect to the parent
- `ω`: [rad/s] Angular velocity of reference frame with respect to the parent
- `a`: [m/s²] Acceleration of reference frame with respect to the parent
- `α`: [rad/s²] Angular acceleration of reference frame with respect to the parent
- `g`: [m/s²] Constant gravitational field vector
"""
mutable struct Cartesian2D
    parent::Union{Nothing, Cartesian2D}
    children::Vector{Cartesian2D}
    translation::Union{Nothing, Vector{Float64}}
    rotation::Union{Nothing, Vector{Float64}}
    v::Union{Nothing, Vector{Float64}}
    ω::Union{Nothing, Vector{Float64}}
    a::Union{Nothing, Vector{Float64}}
    α::Union{Nothing, Vector{Float64}}
    g::UniformField

    function Cartesian2D(parent::Union{Nothing, Cartesian2D} = nothing,
        translation::Union{Nothing, Vector{Float64}} = [0.0,0.0],
        rotation::Union{Nothing, Vector{Float64}}    = [0.0,0.0], 
        v::Union{Nothing, Vector{Float64}}           = [0.0,0.0],
        ω::Union{Nothing, Vector{Float64}}           = [0.0,0.0],
        a::Union{Nothing, Vector{Float64}}           = [0.0,0.0],
        α::Union{Nothing, Vector{Float64}}           = [0.0,0.0],
        g::UniformField                              = UniformField())

        @assert(all(0 .≤ rotation .≤ 2π), "Elements of rotation matrix must be between 0 and 2π radians")
        
        obj = new(parent, Vector{Cartesian2D}[], translation, rotation, v, ω, a, α, g)

        if isnothing(parent)
            @assert(iszero(translation) && iszero(rotation) && iszero(v) && iszero(ω) && iszero(a) && iszero(α), "World frame cannot be transformed or initialized with properties other than a gravitational field")
            obj.translation = nothing; obj.rotation = nothing; obj.v = nothing; obj.ω = nothing; obj.a = nothing; obj.α = nothing
        else
            if obg.g != UniformField()
                @warn "Parameter `g` will be inherited from parent frame"
            end

            obj.g = parent.g
            push!(parent.children, obj)
        end

        return obj
    end
end