"""
Library to model 1-dimensional, rotational mechanical systems
"""
module Rotational

using ModelingToolkitBase, Symbolics, IfElse, SciCompDSL
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ...Blocks: RealInput, RealOutput
import ...@symcheck

export Flange, Support
include("utils.jl")

export Fixed, Inertia, Spring, Damper, SpringDamper, IdealGear, RotationalFriction
include("components.jl")

export Torque, ConstantTorque, Speed, Position
include("sources.jl")

export AngleSensor, SpeedSensor, TorqueSensor, RelSpeedSensor
include("sensors.jl")

end
