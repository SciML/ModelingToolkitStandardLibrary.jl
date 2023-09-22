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
export Flange, Support
include("../Rotational/utils.jl")

"""
    Fixed(;name, phi0 = 0.0)

Flange fixed in housing at a given angle.

# Connectors:

  - `flange` [Flange](@ref)

# Parameters:

  - `phi0`: [`rad`] Fixed offset angle of housing
"""
@mtkmodel Fixed begin
    @components begin
        flange = Flange()
    end
    @parameters begin
        phi0 = 0.0, [description = "Fixed offset angle of flange"]
    end
    @equations begin
        flange.phi ~ phi0
    end
end
end

export Frame, PartialTwoFrames
include("utils.jl")

export Fixed, Body, FixedTranslation
include("components.jl")

export Revolute
include("joints.jl")
end
