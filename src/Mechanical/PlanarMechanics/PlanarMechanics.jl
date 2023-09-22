"""
Library to model planar mechanical multi-body systems inspired by https://github.com/dzimmer/PlanarMechanics
"""

module PlanarMechanics

using ModelingToolkit, Symbolics, IfElse
using ...Blocks: RealInput, RealOutput
import ...@symcheck

@parameters t
D = Differential(t)

module Rotational
# TODO: figure out how to use Rotational directly
using ModelingToolkit
export Flange, Support, Fixed
include("../Rotational/utils.jl")
end

export Frame, PartialTwoFrames
include("utils.jl")

export Fixed, FixedTranslation
include("components.jl")

end
