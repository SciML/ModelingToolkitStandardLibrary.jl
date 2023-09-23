"""
Library to model planar mechanical multi-body systems inspired by https://github.com/dzimmer/PlanarMechanics
"""

module PlanarMechanics

using ModelingToolkit, Symbolics, IfElse
using ...Blocks: RealInput, RealOutput
import ...@symcheck
import ModelingToolkitStandardLibrary.Mechanical.Rotational

@parameters t
D = Differential(t)

export Frame, PartialTwoFrames
include("utils.jl")

export Fixed, Body, FixedTranslation
include("components.jl")

export Revolute
include("joints.jl")
end
