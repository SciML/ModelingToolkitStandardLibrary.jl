# Linearized Components: LinearizedDiode
# Simplified linear models for small-signal analysis

"""
    LinearizedDiode(; name, Vd = 0.7, Rd = 10.0)

Linearized diode model for simplified circuit analysis.

A piecewise-linear diode model with a fixed forward voltage drop and dynamic
resistance. When forward biased (v > Vd), the diode conducts with resistance Rd.
When reverse biased, the diode blocks current.

# Parameters
- `Vd`: [V] Forward voltage drop (default: 0.7)
- `Rd`: [Ω] Dynamic resistance when conducting (default: 10.0)

# Connectors
- `p` Positive terminal (anode)
- `n` Negative terminal (cathode)

# Equations
```
if v > Vd:
    i = (v - Vd) / Rd    # Forward conducting
else:
    i = 0                 # Reverse blocking
```

# Behavior
- Forward bias (v > Vd): Diode conducts, i = (v - Vd) / Rd
- Reverse bias (v ≤ Vd): Diode blocks, i = 0

# Example
```julia
@named diode = LinearizedDiode(Vd = 0.7, Rd = 10.0)
```
"""
@component function LinearizedDiode(; name, Vd = 0.7, Rd = 10.0)
    @named oneport = OnePort()
    @unpack v, i = oneport

    pars = @parameters begin
        Vd = Vd, [description = "Forward voltage drop [V]"]
        Rd = Rd, [description = "Dynamic resistance [Ω]"]
    end

    # Piecewise-linear behavior:
    # Forward bias: i = (v - Vd) / Rd
    # Reverse bias: i = 0
    eqs = Equation[
        i ~ IfElse.ifelse(v > Vd, (v - Vd) / Rd, 0.0),
    ]

    sys = System(eqs, t, [], pars; name, systems = [])
    return extend(sys, oneport)
end
