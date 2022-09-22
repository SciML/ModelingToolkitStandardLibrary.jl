module MultiBody2D

using ModelingToolkit, Symbolics, IfElse

@parameters t
D = Differential(t)

export Port
include("utils.jl")

export Link
include("components.jl")

end
