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
export Sin, Cos, Tan, Asin, Acos, Atan, Atan2, Sinh, Cosh, Tanh, Exp, Log, Log10
include("math.jl")

export Constant, SinSource, ClockSource, RampSource, StepSource
include("sources.jl")

export Limiter, DeadZone, SlewRateLimiter
include("nonlinear.jl")

export Integrator, Derivative, FirstOrder, SecondOrder, StateSpace
export PI, LimPI, PID
include("continuous.jl")

end