"""
Library to model planar mechanical multi-body systems inspired by https://github.com/dzimmer/PlanarMechanics
"""

module PlanarMechanics

import ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit, Symbolics, IfElse
using ...Blocks: RealInput, RealOutput
import ...@symcheck

export Frame, FrameResolve, PartialTwoFrames, ZeroPosition
include("utils.jl")

export Fixed, Body, FixedTranslation, Spring, Damper, SpringDamper
include("components.jl")

export Revolute, Prismatic
include("joints.jl")

export AbsolutePosition,
       RelativePosition, AbsoluteVelocity, RelativeVelocity, AbsoluteAcceleration,
       RelativeAcceleration, connect_sensor
include("sensors.jl")
end
