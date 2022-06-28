"""
    Ground(; name)

Ground node with the potential of zero and connector `g`. Every circuit must have one ground
node.

# Connectors:
- `g`
"""
function Ground(; name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], []; systems = [g], name = name)
end

"""
    Resistor(; name, R)

Creates an ideal Resistor following Ohm's Law.

# States:
See [OnePort](@ref)

# Connectors:
- `p` Positive pin
- `n` Negative pin

# Parameters:
- `R`: [`Ω`] Resistance
"""
function Resistor(; name, R)
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters R = R
    eqs = [
        v ~ i * R,
    ]
    extend(ODESystem(eqs, t, [], pars; name = name), oneport)
end

"""
    Conductor(;name, G)

Creates an ideal conductor.

# States:
See [OnePort](@ref)

# Connectors:
- `p` Positive pin
- `n` Negative pin

# Parameters:
- `G`: [`S`] Conductance
"""
function Conductor(; name, G)
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters G = G
    eqs = [
        i ~ v * G,
    ]
    extend(ODESystem(eqs, t, [], pars; name = name), oneport)
end

"""
    Capacitor(; name, C)


Creates an ideal capacitor.

# States:
- `v(t)`: [`V`]
  The voltage across the capacitor, given by `D(v) ~ p.i / C`

# Connectors:
- `p` Positive pin
- `n` Negative pin

# Parameters:
- `C`: [`F`] Capacitance
- `v_start`: [`V`] Initial voltage of capacitor
"""
function Capacitor(; name, C, v_start = 0.0)
    @named oneport = OnePort(; v_start = v_start)
    @unpack v, i = oneport
    pars = @parameters C = C
    eqs = [
        D(v) ~ i / C,
    ]
    extend(ODESystem(eqs, t, [], pars; name = name), oneport)
end

"""
    Inductor(; name, L)

Creates an ideal Inductor.

# States:
See [OnePort](@ref)

# Connectors:
- `p` Positive pin
- `n` Negative pin

# Parameters:
- `L`: [`H`] Inductance
- `i_start`: [`A`] Initial current through inductor
"""
function Inductor(; name, L, i_start = 0.0)
    @named oneport = OnePort(; i_start = i_start)
    @unpack v, i = oneport
    pars = @parameters L = L
    eqs = [
        D(i) ~ 1 / L * v,
    ]
    extend(ODESystem(eqs, t, [], pars; name = name), oneport)
end

"""
    IdealOpAmp(; name)

Ideal operational amplifier (norator-nullator pair).
The ideal OpAmp is a two-port. The left port is fixed to `v1 = 0` and `i1 = 0` (nullator).
At the right port both any voltage `v2` and any current `i2` are possible (norator).

# States:
See [TwoPort](@ref)

# Connectors:
- `p1` Positive pin (left port)
- `p2` Positive pin (right port)
- `n1` Negative pin (left port)
- `n2` Negative pin (right port)
"""
function IdealOpAmp(; name)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    eqs = [v1 ~ 0
           i1 ~ 0]
    extend(ODESystem(eqs, t, [], [], name = name), twoport)
end

"""
    Short(; name)

Short is a simple short cut branch. That means the voltage drop between both pins is zero.

# States:
See [OnePort](@ref)

# Connectors:
- `p` Positive pin
- `n` Negative pin
"""
function Short(; name)
    @named oneport = OnePort()
    @unpack v, i = oneport
    eqs = [v ~ 0]
    extend(ODESystem(eqs, t, [], []; name = name), oneport)
end
