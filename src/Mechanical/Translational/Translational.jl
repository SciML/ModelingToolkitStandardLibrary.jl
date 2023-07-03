"""
Library to model 1-dimensional, translational mechanical systems
"""
module Translational

using ModelingToolkit, Symbolics
using ModelingToolkit: getdefault

using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using IfElse: ifelse

@parameters t
D = Differential(t)

export MechanicalPort
include("utils.jl")

export Mass, Spring, Damper, Fixed
include("components.jl")

export Force, Position
include("sources.jl")

export ForceSensor, PositionSensor
include("sensors.jl")

end
