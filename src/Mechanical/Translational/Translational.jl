"""
Library to model 1-dimensional, translational mechanical systems
"""
module Translational

using ModelingToolkit, Symbolics
using ModelingToolkit: getdefault

using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using IfElse: ifelse
using ...DynamicQuantities: @u_str

@parameters t [unit = u"s"]
D = Differential(t)

export MechanicalPort
include("utils.jl")

export Mass, Spring, Damper, Fixed
include("components.jl")

export Force, Position, Velocity, Acceleration
include("sources.jl")

export ForceSensor, PositionSensor
include("sensors.jl")

end
