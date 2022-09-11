"""
Library to model 1-dimensional, translational mechanical systems
"""
module Translational

using ModelingToolkit, Symbolics, IfElse
using ...Blocks: RealInput, RealOutput

@parameters t
D = Differential(t)

export Port
include("utils.jl")

export Body, Spring
export SpringDamperBoundary, ForcedBoundary
include("components.jl")

#export Torque
#include("sources.jl")

end
