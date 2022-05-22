"""
```julia
function Ground(; name)
```

Ground node with the potential of zero and connector `g`. Every circuit must have one ground
node.

# Connectors
- `g`
"""
function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], []; systems=[g], name=name)
end

"""
```julia
function Resistor(; name, R)
```

Creates an ideal Resistor following Ohm's Law.

# States
- `v(t)`: [`V`]
  The voltage across the resistor, given by `p.i * R`

# Connectors
- `p`
  Positive pin
- `n`
  Negative pin

# Parameters: 
- `R`: [`Ω`]
  Resistance
"""
function Resistor(;name, R)
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters R=R
    eqs = [
        v ~ i * R
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    Conductor(;name, G)

Ideal linear electrical conductor.

# States
- see [`OnePort`](@ref)

# Connectors
- `p` Positive pin
- `n` Negative pin

# Parameters: 
- `G`: [`S`] Conductance
"""
function Conductor(;name, G)
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters G=G
    eqs = [
        i ~ v * G
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
```julia
function Capacitor(; name, C)
```

Creates an ideal Capacitor.

# States
- `v(t)`: [`V`]
  The voltage across the capacitor, given by `D(v) ~ p.i / C`

# Connectors
- `p`
  Positive pin
- `n`
  Negative pin

# Parameters:
- `C`: [`F`]
  Capacitance
- `v_start`: [`V`]
  Initial voltage of capacitor
"""
function Capacitor(;name, C, v_start=0.0) 
    @named oneport = OnePort(;v_start=v_start)
    @unpack v, i = oneport
    pars = @parameters C=C
    eqs = [
        D(v) ~ i / C
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
```julia
function Inductor(; name, L)
```

Creates an ideal Inductor.

# States
- `v(t)`: [`V`]
  The voltage across the inductor, given by `D(p.i) ~ v / L`

# Connectors
- `p`
  Positive pin
- `n`
  Negative pin

# Parameters:
- `L`: [`H`]
  Inductance
- `i_start`: [`A`]
  Initial current through inductor
"""
function Inductor(;name, L, i_start=0.0)
    @named oneport = OnePort(;i_start=i_start)
    @unpack v, i = oneport
    pars = @parameters L=L
    eqs = [
        D(i) ~ 1 / L * v
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
```julia
function IdealOpAmp(; name)
```

Ideal operational amplifier (norator-nullator pair).
The ideal OpAmp is a two-port. The left port is fixed to `v1 = 0` and `i1 = 0` (nullator). 
At the right port both any voltage `v2` and any current `i2` are possible (norator).

# States
- `v1(t)`: [`V`]
  Voltage of left port
- `v2(t)`: [`V`]
  Voltage of right port
- `i1(t)`: [`A`]
  Current of left port
- `i2(t)`: [`A`]
  Current of right port

# Connectors
- `p1`
  Positive pin (left port)
- `p2`
  Positive pin (right port)
- `n1`
  Negative pin (left port)
- `n2`
  Negative pin (right port)
"""
function IdealOpAmp(;name)
    @named p1 = Pin()
    @named p2 = Pin()
    @named n1 = Pin()
    @named n2 = Pin()
    sts = @variables begin
        v1(t) 
        v2(t) 
        i1(t) 
        i2(t)
    end
    eqs = [
        v1 ~ p1.v - n1.v
        v2 ~ p2.v - n2.v
        0 ~ p1.i + n1.i
        0 ~ p2.i + n2.i
        i1 ~ p1.i
        i2 ~ p2.i
        v1 ~ 0
        i1 ~ 0
    ]
    ODESystem(eqs, t, sts, [], systems=[p1, p2, n1, n2], name=name)
end
