"""
Library to model 1-dimensional, rotational mechanical systems
"""
module Rotational

using ModelingToolkit, Symbolics, IfElse
using ...Blocks: RealInput, RealOutput
import ..@symcheck

@parameters t
D = Differential(t)

export Flange, Support
include("utils.jl")

export Fixed, Inertia, Spring, Damper, IdealGear, RotationalFriction
include("components.jl")

export Torque, ConstantTorque, Speed
include("sources.jl")

export AngleSensor, SpeedSensor, TorqueSensor, RelSpeedSensor
include("sensors.jl")

end
