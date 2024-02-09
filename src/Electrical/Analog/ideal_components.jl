"""
    Ground(; name)

Ground node with the potential of zero and connector `g`. Every circuit must have one ground
node.

# Connectors:

  - `g`
"""
@mtkmodel Ground begin
    @components begin
        g = Pin()
    end
    @equations begin
        g.v ~ 0
    end
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

  - `R`: [`Ohm`] Resistance
"""
@mtkmodel Resistor begin
    @extend v, i = oneport = OnePort()
    @parameters begin
        R, [description = "Resistance", unit = u"Ω"]
    end
    @equations begin
        v ~ i * R
    end
end

"""
    Conductor(; name, G)

Creates an ideal conductor.

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin

# Parameters:

  - `G`: [`S`] Conductance
"""
@mtkmodel Conductor begin
    @extend v, i = oneport = OnePort()
    @parameters begin
        G, [description = "Conductance", unit = u"S"]
    end
    @equations begin
        i ~ v * G
    end
end

"""
    Capacitor(; name, C, v)

Creates an ideal capacitor.
Initial voltage of capacitor can be set with `v` ([`V`])

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin

# Parameters:

  - `C`: [`F`] Capacitance
"""
@mtkmodel Capacitor begin
    @extend v, i = oneport = OnePort(; v = 0.0)
    @parameters begin
        C, [description = "Capacitance",unit = u"F"]
    end
    @equations begin
        D(v) ~ i / C
    end
end

"""
    Inductor(; name, L, i)

Creates an ideal Inductor.
Initial current through inductor can be set with `i` ([`A`]).

# States:

See [OnePort](@ref)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin

# Parameters:

  - `L`: [`H`] Inductance
"""
@mtkmodel Inductor begin
    @extend v, i = oneport = OnePort(; i = 0.0)
    @parameters begin
        L, [description = "Inductance", unit = u"H"]
    end
    @equations begin
        D(i) ~ 1 / L * v
    end
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
@mtkmodel IdealOpAmp begin
    @extend v1, v2, i1, i2 = twoport = TwoPort()
    @equations begin
        v1 ~ 0
        i1 ~ 0
    end
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
@mtkmodel Short begin
    @extend v, i = oneport = OnePort()
    @equations begin
        v ~ 0
    end
end

"""
    HeatingResistor(; name, R_ref = 1.0, T_ref = 300.15, alpha = 0)

Temperature dependent electrical resistor

# States

  - See [OnePort](@ref)
  - `R(t)`: [`Ohm`] Temperature dependent resistance `R ~ R_ref*(1 + alpha*(heat_port.T(t) - T_ref))`

# Connectors

  - `p` Positive pin
  - `n` Negative pin

# Parameters:

  - `R_ref`: [`Ω`] Reference resistance
  - `T_ref`: [K] Reference temperature
  - `alpha`: [K⁻¹] Temperature coefficient of resistance
"""
@mtkmodel HeatingResistor begin
    @extend v, i = oneport = OnePort()
    @components begin
        heat_port = HeatPort()
    end
    @parameters begin
        R_ref = 1.0, [description = "Reference resistance", unit = u"Ω"]
        T_ref = 300.15, [description = "Reference temperature", unit = u"K"]
        alpha = 0, [description = "Temperature coefficient of resistance", unit = u"1/K"]
    end
    @variables begin
        R(t) = R_ref
    end
    @equations begin
        R ~ R_ref * (1 + alpha * (heat_port.T - T_ref))
        heat_port.Q_flow ~ -v * i # -LossPower
        v ~ i * R
    end
end
"""
    EMF(; name, k)

Electromotoric force (electric/mechanic transformer)

# States

  - `v(t)`: [`V`] The voltage across component `p.v - n.v`
  - `i(t)`: [`A`] The current passing through positive pin
  - `phi`: [`rad`] Rotation angle (=flange.phi - support.phi)
  - `w`: [`rad/s`] Angular velocity (= der(phi))

# Connectors

  - `p` [Pin](@ref) Positive pin
  - `n` [Pin](@ref) Negative pin
  - `flange` [Flange](@ref) Shaft of EMF shaft
  - `support` [Support](@ref) Support/housing of emf shaft

# Parameters:

  - `k`: [`N⋅m/A`] Transformation coefficient
"""
@mtkmodel EMF begin
    @extend OnePort()
    @components begin
        flange = Flange()
        support = Support()
    end
    @parameters begin
        k, [description = "Transformation coefficient", unit = u"N*m/A"]
    end
    @variables begin
        phi(t) = 0.0, [description = "Rotation angle", unit = u"rad"]
        w(t) = 0.0, [description = "Angular velocity", unit = u"rad/s"]
    end
    @equations begin
        phi ~ flange.phi - support.phi
        D(phi) ~ w
        k * w ~ v
        flange.tau ~ -k * i
    end
end
