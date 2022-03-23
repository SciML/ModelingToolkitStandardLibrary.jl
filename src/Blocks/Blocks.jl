"""
The module `Blocks` contains common input-output components, referred to as blocks.

In general, input-output blocks follow the convention 
```
     ┌───────────┐
 u   │  ẋ=f(x,u) │  y
────►│  y=g(x,u) ├────►
     │           │
     └───────────┘
```
where `u` are inputs, `x` are state variables and `y` are outputs. `x,u,y` are all implemented as `@variables` internally, `u` are marked as `[input=true]` and `y` are marked `[output=true]`.
"""
module Blocks
using ModelingToolkit, Symbolics, IfElse, OrdinaryDiffEq

@parameters t
D = Differential(t)

include("utils.jl")

export Gain, Sum
include("math.jl")

export Saturation, DeadZone
include("nonlinear.jl")

export Constant, Integrator, Derivative, FirstOrder, SecondOrder, PID, StateSpace
include("continuous.jl")

end