"""
Library to model 1-dimensional, translational mechanical components.
"""
module TranslationalModelica

using ModelingToolkit, Symbolics, IfElse
using ...Blocks: RealInput, RealOutput
using ...DynamicQuantities: @u_str

@parameters t [unit = u"s"]
D = Differential(t)

export Flange
include("utils.jl")

export Fixed, Mass, Spring, Damper
include("components.jl")

export Force
include("sources.jl")

end
