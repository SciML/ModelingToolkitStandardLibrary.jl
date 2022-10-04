module MultiBody2D

using ModelingToolkit, Symbolics, IfElse
using ..Translational

@parameters t
D = Differential(t)

export Link
include("components.jl")

end
