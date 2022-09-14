"""
Library to model 1-dimensional, translational mechanical systems
"""
module Translational

using ModelingToolkit, Symbolics, IfElse

@parameters t
D = Differential(t)

export Port
include("utils.jl")

export Body, Spring, Damper
include("components.jl")

#export Torque
#include("sources.jl")

end
