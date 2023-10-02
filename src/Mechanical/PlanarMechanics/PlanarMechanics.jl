"""
Library to model planar mechanical multi-body systems inspired by https://github.com/dzimmer/PlanarMechanics
"""

module PlanarMechanics

import ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica
using ModelingToolkit, Symbolics, IfElse
using ...Blocks: RealInput, RealOutput
import ...@symcheck

@parameters t
D = Differential(t)

export Frame, FrameResolve, PartialTwoFrames, ZeroPosition
include("utils.jl")

export Fixed, Body, FixedTranslation
include("components.jl")

export Revolute, Prismatic
include("joints.jl")

export AbsolutePosition,
    RelativePosition, AbsoluteVelocity, RelativeVelocity, AbsoluteAcceleration,
    RelativeAcceleration
include("sensors.jl")
end
