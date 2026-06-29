"""
Library to model 1-dimensional, translational mechanical systems
"""
module Translational

using ModelingToolkitBase: ModelingToolkitBase, @component, @connector, @named,
    @parameters, Flow, System, compose, t_nounits as t, D_nounits as D
using Symbolics: Symbolics, @variables, Equation

using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput

export MechanicalPort
include("utils.jl")

export Mass, Spring, Damper, Fixed
include("components.jl")

export Force, Position, Velocity, Acceleration
include("sources.jl")

export ForceSensor, PositionSensor, AccelerationSensor
include("sensors.jl")

end
