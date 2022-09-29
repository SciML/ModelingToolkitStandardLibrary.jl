"""
The module `Blocks` contains common input-output components, referred to as blocks.
"""
module Blocks
using ModelingToolkit, Symbolics
using IfElse: ifelse

@parameters t
D = Differential(t)

export RealInput, RealOutput, SISO
include("utils.jl")

export Gain, Sum, MatrixGain, Feedback, Add, Add3, Product, Division
export Abs, Sign, Sqrt, Sin, Cos, Tan, Asin, Acos, Atan, Atan2, Sinh, Cosh, Tanh, Exp
export Log, Log10
include("math.jl")

export Constant, Sine, Cosine, ContinuousClock, Ramp, Step, ExpSine, Square, Triangular
include("sources.jl")

export Limiter, DeadZone, SlewRateLimiter
include("nonlinear.jl")

export Integrator, Derivative, FirstOrder, SecondOrder, StateSpace
export PI, LimPI, PID, LimPID
include("continuous.jl")

export AnalysisPoint, get_sensitivity, get_comp_sensitivity,
       get_looptransfer, open_loop
include("analysis_points.jl")

end
