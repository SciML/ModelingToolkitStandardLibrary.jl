"""
Library to model 1-dimensional, translational mechanical components.
"""
module TranslationalPosition

using ModelingToolkit, Symbolics, IfElse
using ...Blocks: RealInput, RealOutput

@parameters t
D = Differential(t)

export Flange
include("utils.jl")

export Fixed, Mass, Spring, Damper
include("components.jl")

export Force
include("sources.jl")

end
