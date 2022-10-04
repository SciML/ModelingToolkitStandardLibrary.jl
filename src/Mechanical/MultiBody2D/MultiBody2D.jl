module MultiBody2D

using ModelingToolkit, Symbolics, IfElse
using ..Translational

@parameters t
D = Differential(t)

export Link, RevoluteJoint, MultiBody2Translational
include("components.jl")

export RigidBody2DPort
include("utils.jl")

end
