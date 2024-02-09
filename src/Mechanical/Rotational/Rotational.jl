"""
Library to model 1-dimensional, rotational mechanical systems
"""
module Rotational

using ModelingToolkit, Symbolics, IfElse
using ...Blocks: RealInput, RealOutput
import ...@symcheck
using ...DynamicQuantities: @u_str
import ...Wb
import ...rad

@parameters t [unit = u"s"]
D = Differential(t)

export Flange, Support
include("utils.jl")

export Fixed, Inertia, Spring, Damper, SpringDamper, IdealGear, RotationalFriction
include("components.jl")

export Torque, ConstantTorque, Speed
include("sources.jl")

export AngleSensor, SpeedSensor, TorqueSensor, RelSpeedSensor
include("sensors.jl")


end
