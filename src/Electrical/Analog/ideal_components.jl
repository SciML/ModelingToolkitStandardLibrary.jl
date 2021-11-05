"""
```julia
function Ground(;name)
```

Ground node with the potential of zero and connector `g`

# Connectors
- `g`
"""
function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], [], systems=[g], name=name)
end

"""
```julia
function Resistor(;name, R = 1.0)
```

Creates an ideal Resistor following Ohm's Law.

# Observables
- `R`
  Resistance (negative, zero, positive)

# States
- `v(t)`
  The voltage across the resistor, given by `p.i * R`

# Connectors
- `p`
  Positive pin
- `n`
  Negative pin
"""
function Resistor(;name, R = 1.0)
    val = R
    
    @named p = Pin()
    @named n = Pin()
    @parameters R
    @variables v(t)

    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           v ~ p.i * R
          ]
    ODESystem(eqs, t, [v], [R], systems=[p, n], defaults=Dict(R => val), name=name)
end

"""
```julia
function Capacitor(; name, C = 1.0)
```

Creates an ideal Capacitor.

# Observables
- `C`
  Capacitance (zero or positive)

# States
- `v(t)`
  The voltage across the capacitor, given by `D(v) ~ p.i / C`

# Connectors
- `p`
  Positive pin
- `n`
  Negative pin
"""
function Capacitor(; name, C = 1.0)
    val = C

    @named p = Pin()
    @named n = Pin()
    @parameters C
    @variables v(t)

    D = Differential(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           D(v) ~ p.i / C
          ]
    ODESystem(eqs, t, [v], [C], systems=[p, n], defaults=Dict(C => val), name=name)
end

"""
```julia
function Inductor(; name, L = 1.0)
```

Creates an ideal Inductor.

# Observables
- `L`
  Inductance (zero or positive)

# States
- `v(t)`
  The voltage across the inductor, given by `D(p.i) ~ v / L`

# Connectors
- `p`
  Positive pin
- `n`
  Negative pin
"""
function Inductor(; name, L = 1.0)
    val = L

    @named p = Pin()
    @named n = Pin()
    @parameters L
    @variables v(t)

    D = Differential(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           D(p.i) ~ v / L
          ]
    ODESystem(eqs, t, [v], [L], systems=[p, n], defaults=Dict(L => val), name=name)
end

"""
```julia
function IdealOpAmp(; name)
```

Creates an ideal Operational Amplifier.

# States
- `v1(t)`
  Voltage of left port
- `v2(t)`
  Voltage of right port
- `i1(t)`
  Current of left port
- `i2(t)`
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
function IdealOpAmp(; name)
    @named p1 = Pin()
    @named p2 = Pin()
    @named n1 = Pin()
    @named n2 = Pin()
    @variables v1(t) v2(t) # u"v"
    @variables i1(t) i2(t) # u"A"

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
    ODESystem(eqs, t, [i1, i2, v1, v2], [], systems=[p1, p2, n1, n2], name=name)
end
