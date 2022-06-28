"""
Library to model 1-dimensional, rotational mechanical systems
"""
module Rotational

using ModelingToolkit, Symbolics
using ...Blocks: RealInput, RealOutput

@parameters t
D = Differential(t)

export Flange, Support
include("utils.jl")

export Fixed, Inertia, Spring, Damper, IdealGear, RotationalFriction
include("components.jl")

export Torque, Speed
include("sources.jl")

export AngleSensor, SpeedSensor, TorqueSensor, RelSpeedSensor
include("sensors.jl")

end
