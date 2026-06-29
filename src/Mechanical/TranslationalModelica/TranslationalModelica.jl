"""
Library to model 1-dimensional, translational mechanical components.
"""
module TranslationalModelica

using IfElse: IfElse
using ModelingToolkitBase: ModelingToolkitBase, @component, @connector, @named,
    @parameters, @unpack, Flow, System, compose, extend, t_nounits as t, D_nounits as D
using Symbolics: Symbolics, @variables, Equation
using ...Blocks: RealInput

export Flange
include("utils.jl")

export Fixed, Mass, Spring, Damper, SpringDamper
include("components.jl")

export Force, Position
include("sources.jl")

end
