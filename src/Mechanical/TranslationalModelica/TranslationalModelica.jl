"""
Library to model 1-dimensional, translational mechanical components.
"""
module TranslationalModelica

using ModelingToolkitBase, Symbolics, IfElse, SciCompDSL
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ...Blocks: RealInput, RealOutput

export Flange
include("utils.jl")

export Fixed, Mass, Spring, Damper, SpringDamper
include("components.jl")

export Force, Position
include("sources.jl")

end
