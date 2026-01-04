"""
    Ground(; name)

Ground node with the potential of zero and connector `g`. Every circuit must have one ground
node.

# Connectors:

  - `g`
"""
@component function Ground(; name)
    pars = @parameters begin
    end

    systems = @named begin
        g = Pin()
    end

    vars = @variables begin
    end

    equations = Equation[
        g.v ~ 0,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Resistor(; name, R = 1.0, T_ref = 300.15, alpha = 0, T_dep = false)

Generic resistor with optional temperature dependency.

# States:

  - See [OnePort](@ref)
  - `R(t)`: [`Ω`] Resistance (temperature dependent if `T_dep = true`)

# Connectors:

  - `p` Positive pin
  - `n` Negative pin
  - `heat_port` [HeatPort](@ref) (only if `T_dep = true`) Heat port to model the temperature dependency

# Parameters:

  - `R`: [`Ω`] Reference resistance
  - `T_ref`: [K] Reference temperature
  - `alpha`: [K⁻¹] Temperature coefficient of resistance
  - `T_dep`: [bool] Temperature dependency
"""
@component function Resistor(; T_dep = false, R = 1.0, T_ref = 300.15, alpha = 0.0, name)
    @named oneport = OnePort()
    @unpack v, i = oneport

    pars = @parameters begin
        R = R, [description = "Reference resistance"]
        T_ref = T_ref, [description = "Reference temperature"]
        alpha = alpha, [description = "Temperature coefficient of resistance"]
    end

    systems = if T_dep
        @named begin
            heat_port = HeatPort()
        end
    else
        @named begin
        end
    end

    vars = if T_dep
        @variables begin
            R_T(t), [description = "Temperature-dependent resistance"]
        end
    else
        @variables begin
        end
    end

    equations = if T_dep
        Equation[
            R_T ~ R * (1 + alpha * (heat_port.T - T_ref)),
            heat_port.Q_flow ~ -v * i,
            v ~ i * R_T,
        ]
    else
        Equation[
            v ~ i * R,
        ]
    end

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, oneport)
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
@component function Conductor(; G = nothing, name)
    @named oneport = OnePort()
    @unpack v, i = oneport

    pars = @parameters begin
        G = G, [description = "Conductance"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        i ~ v * G,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, oneport)
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
@component function Capacitor(; C = nothing, v = nothing, name)
    @named oneport = OnePort(; v)
    @unpack v, i = oneport

    pars = @parameters begin
        C = C, [description = "Capacitance"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        D(v) ~ i / C,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, oneport)
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
@component function Inductor(; L = nothing, i = nothing, name)
    @named oneport = OnePort(; i)
    @unpack v, i = oneport

    pars = @parameters begin
        L = L, [description = "Inductance"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        D(i) ~ 1 / L * v,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, oneport)
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
@component function IdealOpAmp(; name)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    pars = @parameters begin
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        v1 ~ 0,
        i1 ~ 0,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, twoport)
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
@component function Short(; name)
    @named oneport = OnePort()
    @unpack v, i = oneport

    pars = @parameters begin
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        v ~ 0,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, oneport)
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
@component function EMF(; k = nothing, name)
    @named oneport = OnePort()
    @unpack v, i = oneport

    pars = @parameters begin
        k = k, [description = "Transformation coefficient"]
    end

    systems = @named begin
        flange = Flange()
        support = Support()
    end

    vars = @variables begin
        phi(t), [guess = 0.0, description = "Rotation Angle"]
        w(t), [guess = 0.0]
    end

    equations = Equation[
        phi ~ flange.phi - support.phi,
        D(phi) ~ w,
        k * w ~ v,
        flange.tau ~ -k * i,
        support.tau ~ -flange.tau,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, oneport)
end

"""
    Diode(; name, Is = 1e-6, n = 1, T = 300.15, T_dep = false)

Generic diode with optional temperature dependency.

# States

    - See [OnePort](@ref)

# Connectors

    - `p` Positive pin
    - `n` Negative pin
    - `port` [HeatPort](@ref) (only if `T_dep = true`) Heat port to model variable temperature dependency

# Parameters:
    
    - `Is`: [`A`] Saturation current
    - `n`: Ideality factor
    - `T`: [K] Constant ambient temperature - only used if T_dep=false
    - `T_dep`: [bool] Temperature dependency
"""
@component function Diode(; T_dep = false, Is = 1.0e-6, n = 1, T = 300.15, v = 0.0, name)
    consts = @constants begin
        k = 1.380649e-23 # Boltzmann constant (J/K)
        q = 1.602176634e-19 # Elementary charge (C)
    end

    @named oneport = OnePort(; v)
    @unpack v, i = oneport

    pars = @parameters begin
        Is = Is, [description = "Saturation current (A)"]
        n = n, [description = "Ideality factor"]
        T = T, [description = "Ambient temperature"]
    end
    pars = [pars; consts]

    systems = if T_dep
        @named begin
            port = HeatPort()
        end
    else
        @named begin
        end
    end

    vars = if T_dep
        @variables begin
            Vt(t), [description = "Thermal voltage"]
        end
    else
        @variables begin
        end
    end

    equations = if T_dep
        Equation[
            Vt ~ k * port.T / q,  # Thermal voltage equation
            i ~ Is * (exp(v / (n * Vt)) - 1),  # Shockley diode equation with temperature dependence
            port.Q_flow ~ -v * i,  # -LossPower
        ]
    else
        Equation[
            i ~ Is * (exp(v * q / (n * k * T)) - 1),  # Shockley diode equation
        ]
    end

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, oneport)
end

"""
    VariableResistor(; name, R_ref = 1.0, T_ref = 300.15, R_const = 1e-3, T_dep = false)

Variable resistor with optional temperature dependency.

The total resistance R ∈ [R_const, R_const + R_ref], where pos is the
position of the wiper and R_ref is the variable resistance between p and n.
The total resistance is then:

R = R_const + pos * R_ref

If T_dep is true, then R also depends on the temperature of the heat port with
temperature coefficient alpha. The total resistance is then:

R = R_const + pos * R_ref * (1 + alpha * (port.T - T_ref))

# States

    - See [OnePort](@ref)
    - `pos(t)`: Position of the wiper (normally 0-1)
    - `R(t)`: Resistance

# Connectors

        - `p` Positive pin
        - `n` Negative pin
        - `position` RealInput to set the position of the wiper
        - `port` [HeatPort](@ref) Heat port to model the temperature dependency

# Parameters

        - `R_ref`: [`Ω`] Resistance at temperature T_ref when fully closed (pos=1.0)
        - `T_ref`: [K] Reference temperature
        - `R_const`: [`Ω`] Constant resistance between p and n
        - `T_dep`: Temperature dependency
        - `alpha`: [K⁻¹] Temperature coefficient of resistance
        - `enforce_bounds`: Enforce bounds for the position of the wiper (0-1)
"""
@component function VariableResistor(; T_dep = false, enforce_bounds = true, R_ref = 1.0, T_ref = 300.15, R_const = 1.0e-3, alpha = 1.0e-3, name)
    @named oneport = OnePort()
    @unpack v, i = oneport

    pars = if T_dep
        @parameters begin
            R_ref = R_ref, [description = "Resistance at temperature T_ref when fully closed (pos=1.0) (Ω)"]
            T_ref = T_ref, [description = "Reference temperature (K)"]
            R_const = R_const, [description = "Constant resistance between p and n (Ω)"]
            alpha = alpha, [description = "Temperature coefficient of resistance (K^-1)"]
        end
    else
        @parameters begin
            R_ref = R_ref, [description = "Resistance at temperature T_ref when fully closed (pos=1.0) (Ω)"]
            T_ref = T_ref, [description = "Reference temperature (K)"]
            R_const = R_const, [description = "Constant resistance between p and n (Ω)"]
        end
    end

    systems = if T_dep
        @named begin
            position = RealInput()
            port = HeatPort()
        end
    else
        @named begin
            position = RealInput()
        end
    end

    vars = @variables begin
        pos(t), [description = "Position of the wiper (normally 0-1)"]
        R(t), [description = "Resistance (Ω)"]
    end

    conditional_eqs = if T_dep
        Equation[
            port.Q_flow ~ -v * i,  # -LossPower
            R ~ R_const + pos * R_ref * (1 + alpha * (port.T - T_ref)),
        ]
    else
        Equation[
            R ~ R_const + pos * R_ref,
        ]
    end

    equations = [
        conditional_eqs...
        pos ~ (enforce_bounds ? clamp(position.u, 0, 1) : position.u)
        v ~ i * R
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, oneport)
end
