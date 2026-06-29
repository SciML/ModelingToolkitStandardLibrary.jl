"""
Library to model 1-dimensional, rotational mechanical systems
"""
module Rotational

using IfElse: IfElse
using ModelingToolkitBase: ModelingToolkitBase, @component, @connector, @named,
    @parameters, @unpack, Flow, System, compose, extend, t_nounits as t, D_nounits as D
using Symbolics: Symbolics, @variables, Equation
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
