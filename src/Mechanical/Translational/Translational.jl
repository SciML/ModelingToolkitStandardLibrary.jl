module Translational

using ModelingToolkit, Symbolics, IfElse, OrdinaryDiffEq
using OffsetArrays
using ...Blocks: RealInput, RealOutput

@parameters t
D = Differential(t)

export Flange
include("utils.jl")

export Fixed, Inertia, Spring, Damper, IdealGear
include("components.jl")

export Force
include("sources.jl")

end