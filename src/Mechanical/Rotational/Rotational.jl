"""
Library to model 1-dimensional, rotational mechanical systems
"""
module Rotational

using ModelingToolkit, Symbolics, IfElse
using ...Blocks: RealInput, RealOutput

@parameters t
D = Differential(t)

export Flange
include("utils.jl")

export Fixed, Inertia, Spring, Damper, IdealGear
include("components.jl")

export Torque
include("sources.jl")

end
