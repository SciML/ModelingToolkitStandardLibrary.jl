"""
Library to model 1-dimensional, rotational mechanical systems
"""
module Rotational

using ModelingToolkit, Symbolics, IfElse, OrdinaryDiffEq
using OffsetArrays
using ...Blocks: RealInput, RealOutput

@parameters t
D = Differential(t)

export Flange, Support
include("utils.jl")

export Fixed, Inertia, Spring, Damper, IdealGear, RotationalFriction
include("components.jl")

export Torque
include("sources.jl")

end