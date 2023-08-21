module MultiBody2D

using ModelingToolkit, Symbolics, IfElse
using ..TranslationalPosition

@parameters t
D = Differential(t)

export Link
include("components.jl")

end
