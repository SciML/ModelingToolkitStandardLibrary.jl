const ie = IfElse.ifelse

"""
    Saturation(; y_max, y_min=-y_max, name)

The `Saturation` IO block limits the output between `y_min` and `y_max`, equivalent to
`y ~ clamp(u, y_min, y_max)`.

Keywords: limiter, sat, actuator model
"""
function Saturation(; y_max, y_min=y_max > 0 ? -y_max : -Inf, name)
    if !ModelingToolkit.isvariable(y_max)
        y_max ≥ y_min || throw(ArgumentError("y_min must be smaller than y_max"))
    end
    @named u = RealInput()
    @named y = RealOutput()
    pars = @parameters y_max=y_max y_min=y_min
    eqs = [
        # The equation below is equivalent to y ~ clamp(u, y_min, y_max)
        y.u ~ ie(u.u > y_max, y_max, ie( (y_min < u.u) & (u.u < y_max), u.u, y_min))
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [y, u])
end

"""
    DeadZone(; u_max, u_min=-u_max, name)

A dead zone is a band within which the output is zero.
Outside of the dead zone, the output changes linearly starting from zero at the band edge.
```
       y▲
        │     /
        │    /
  u_min │   /
─────|──┼──|───────► u
    /   │   u_max
   /    │
  /     │
```
"""
function DeadZone(; u_max, u_min=-u_max, name)
    if !ModelingToolkit.isvariable(u_max)
        u_max ≥ u_min || throw(ArgumentError("u_min must be smaller than u_max"))
    end
    @named u = RealInput()
    @named y = RealOutput()
    pars = @parameters u_max=u_max u_min=u_min
    eqs = [
        y.u ~ ie(u.u > u_max, u.u - u_max, ie(u.u < u_min, u.u - u_min, 0))
    ]
    compose(ODESystem(eqs, t, [], pars; name=name), [y, u])
end
