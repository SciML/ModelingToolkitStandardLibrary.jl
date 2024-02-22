"""
Library to model 1-dimensional, translational mechanical components.
"""
module TranslationalModelica

using ModelingToolkit, Symbolics, IfElse
using ModelingToolkit: t, D
using DynamicQuantities

using ...Blocks: RealInput, RealOutput

export Flange
include("utils.jl")

export Fixed, Mass, Spring, Damper
include("components.jl")

export Force
include("sources.jl")

end
