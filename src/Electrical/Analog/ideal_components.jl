"""
    Ground(; name)

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
    Resistor(; name, R = 1.0)

Creates an ideal Resistor following Ohm's Law.

# States
See [OnePort](@ref)

# Connectors
- `p` Positive pin
- `n` Negative pin

# Parameters: 
- `R`: [`Ω`] Resistance
"""
function Resistor(;name, R=1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    pars = @parameters R=R
    eqs = [
        v ~ i * R
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    Capacitor(; name, C = 1.0)

Creates an ideal Capacitor.

# States
See [OnePort](@ref)

# Connectors
- `p` Positive pin
- `n` Negative pin

# Parameters:
- `C`: [`F`] Capacitance
- `v_start`: [`V`] Initial voltage of capacitor
"""
function Capacitor(;name, C=1.0, v_start=0.0) 
    @named oneport = OnePort(;v_start=v_start)
    @unpack v, i = oneport
    pars = @parameters C=C
    eqs = [
        D(v) ~ i / C
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    Inductor(; name, L = 1.0)

Creates an ideal Inductor.

# States
See [OnePort](@ref)

# Connectors
- `p` Positive pin
- `n` Negative pin

# Parameters:
- `L`: [`H`] Inductance
- `i_start`: [`A`] Initial current through inductor
"""
function Inductor(;name, L=1.0e-6, i_start=0.0)
    @named oneport = OnePort(;i_start=i_start)
    @unpack v, i = oneport
    pars = @parameters L=L
    eqs = [
        D(i) ~ 1 / L * v
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), oneport)
end

"""
    IdealOpAmp(; name)

Ideal operational amplifier (norator-nullator pair).
The ideal OpAmp is a two-port. The left port is fixed to `v1 = 0` and `i1 = 0` (nullator). 
At the right port both any voltage `v2` and any current `i2` are possible (norator).

# States
- `v1(t)`: [`V`] Voltage of left port
- `v2(t)`: [`V`] Voltage of right port
- `i1(t)`: [`A`] Current of left port
- `i2(t)`: [`A`] Current of right port

# Connectors
- `p1` Positive pin (left port)
- `p2` Positive pin (right port)
- `n1` Negative pin (left port)
- `n2` Negative pin (right port)
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


"""
    HeatingResistor(;name, R_ref=1.0, T_ref=300.15, alpha=0)

Temperature dependent electrical resistor

# States
- See [OnePort](@ref)
- `R(t)`: [`Ω`] Temperature dependent resistance `R ~ R_ref*(1 + alpha*(heat_port.T(t) - T_ref))`

# Connectors
- `p` Positive pin
- `n` Negative pin

# Parameters: 
- `R_ref`: [`Ω`] Reference resistance 
- `T_ref`: [K] Reference temperature
"""
function HeatingResistor(;name, R_ref=1.0, T_ref=300.15, alpha=0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    @names heat_port = HeatPort()
    pars = @parameters begin
      R_ref=R_ref
      T_ref=T_ref
      alpha=alpha
    end
    @variables R
    eqs = [
      R ~ R_ref*(1 + alpha*(heat_port.T - T_ref))
      heat_port.Q_flow ~ -v * i # -LossPower
      v ~ i * R
    ]
    extend(ODESystem(eqs, t, [R], pars; name=name, systems=[heat_port]), oneport)
end

"""
    EMF(;name, k)

Electromotoric force (electric/mechanic transformer)

# States
- `v(t)`: [`V`] The voltage across component `p.v - n.v`
- `i(t)`: [`A`] The current passing through positive pin
- `phi`: [`rad`] Rotation angle (=flange.phi - support.phi)
- `w`: [`rad/s`] Angular velocity (= der(phi))

# Connectors
- `p` Positive pin
- `n` Negative pin
- `flange` Shaft of EMF shaft
- `support` Support/housing of emf shaft

# Parameters: 
- `k`: [`N⋅m/A`] Transformation coefficient 
"""
function EMF(;name, k)
    @named p = Pin()
    @named n = Pin()
    @named flange = Flange()
    @named support = Support()
    @parameters k=k
    @variables v(t) i(t) phi(t) w(t)
    eqs = [
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
        phi ~ flange.phi - support.phi
        D(phi) ~ w
        k*w ~ v
        flange.tau ~ -k*i
    ]
    ODESystem(eqs, t, [v, i, phi, w], [k]; name=name, systems=[p, n, flange, support])
end