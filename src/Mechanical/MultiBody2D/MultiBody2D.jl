module MultiBody2D

using ModelingToolkit, Symbolics, IfElse
using ..Translational
using ...DynamicQuantities: @u_str

@parameters t [unit = u"s"]
D = Differential(t)

export Link
include("components.jl")

end
