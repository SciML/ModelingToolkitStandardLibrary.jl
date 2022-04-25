_clamp(u, u_min, u_max) = max(min(u, u_max), u_min)
_dead_zone(u, u_min, u_max) = ifelse(u > u_max, u - u_max, ifelse(u < u_min, u - u_min, 0))

"""
Limit the range of a signal.

# Parameters:
- `y_max`: Maximum of output signal
- `y_min`: Minimum of output signal
"""
function Limiter(;name, y_max, y_min=y_max > 0 ? -y_max : -Inf)
    y_max ≥ y_min || throw(ArgumentError("`y_min` must be smaller than `y_max`"))
    @named siso = SISO()
    @unpack u, y = siso
    pars = @parameters y_max=y_max y_min=y_min
    eqs = [
        y ~ _clamp(u, y_min, y_max)
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
        y ~ _dead_zone(u, u_min, u_max)
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), siso)
end

"""
    SlewRateLimiter(;name, rising=1, falling=-rising, Td=0.001, y_start=0.0)
    
Limits the slew rate of a signal.

# Parameters:
- `Rising`: Maximum rising slew rate
- `falling`: Maximum falling slew rate
- `Td`: Derivative time constant
"""
function SlewRateLimiter(;name, rising=1, falling=-rising, Td=0.001, y_start=0.0)
    rising ≥ falling || throw(ArgumentError("`rising` must be smaller than `falling`"))
    Td > 0 || throw(ArgumentError("Time constant `Td` must be strictly positive"))
    @named siso = SISO(y_start=y_start)
    @unpack u, y = siso
    pars = @parameters rising=rising falling=falling
    eqs = [
        D(y) ~ max(min((u-y) / Td, rising), falling)
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), siso)
end
