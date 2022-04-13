"""
    Saturation(; y_max, y_min=-y_max, name)

The Limiter block passes its input signal as output signal as long as the input is within the specified upper and lower limits. 
If this is not the case, the corresponding limits are passed as output.
"""
function Saturation(;name, y_max, y_min=y_max > 0 ? -y_max : -Inf)
    y_max ≥ y_min || throw(ArgumentError("`y_min` must be smaller than `y_max`"))
    @named siso = SISO()
    @unpack u, y = siso
    pars = @parameters y_max=y_max y_min=y_min
    eqs = [
        y ~ max(min(u, y_max), y_min)
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), siso)
end

"""
    DeadZone(; u_max, u_min=-u_max, name)

The DeadZone block defines a region of zero output.
If the input is within uMin ... uMax, the output is zero. Outside of this zone, the output is a linear function of the input with a slope of 1.
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
function DeadZone(; name, u_max, u_min=-u_max)
    if !ModelingToolkit.isvariable(u_max)
        u_max ≥ u_min || throw(ArgumentError("`u_min` must be smaller than `u_max`"))
    end
    @named siso = SISO()
    @unpack u, y = siso
    pars = @parameters u_max=u_max u_min=u_min
    eqs = [
        y ~ ifelse(u > u_max, u - u_max, ifelse(u < u_min, u - u_min, 0))
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), siso)
end