"""
The module `Blocks` contains common input-output components, referred to as blocks.
"""
module Blocks
using ModelingToolkit, Symbolics, IfElse, OrdinaryDiffEq
using IfElse: ifelse

@parameters t
D = Differential(t)

export RealInput, RealOutput, SISO
include("utils.jl")

export Gain, Sum, MatrixGain, Sum, Feedback, Add, Product, Division, Abs, Sign, Sqrt
include("math.jl")

export Constant, SinSource
include("sources.jl")

export Saturation, DeadZone
include("nonlinear.jl")

export Integrator, Derivative, FirstOrder, SecondOrder #TODO: , PID, StateSpace
export PI
include("continuous.jl")

end