module MultiBody2D

using ModelingToolkit, Symbolics, IfElse
using ModelingToolkit: t, D

using ..TranslationalPosition

export Link
include("components.jl")

end
