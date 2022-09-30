"""
Library to model 1-dimensional, translational mechanical systems
"""
module Translational

using ModelingToolkit, Symbolics, IfElse

@parameters t
D = Differential(t)

export MechanicalPort
include("utils.jl")

export Body, Spring, Damper, Fixed
include("components.jl")

export PositionSensor
include("sensors.jl")

end
