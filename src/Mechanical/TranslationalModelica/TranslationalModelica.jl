"""
Library to model 1-dimensional, translational mechanical components.
"""
module TranslationalModelica

using ModelingToolkit, Symbolics, IfElse
using ModelingToolkit: t_nounits as t, D_nounits as D
using ...Blocks: RealInput, RealOutput

export Flange, Support
include("utils.jl")

export Fixed, Mass, Spring, Damper
include("components.jl")

export Force
include("sources.jl")

end
