"""
    HeatCapacitor(; name, C, T = 273.15 + 20)

Lumped thermal element storing heat

# States:

  - `T`: [`K`] Temperature of element. It accepts an initial value, which defaults to 273.15 + 20.
  - `der_T`: [`K/s`] Time derivative of temperature

# Connectors:

  - `port`

# Parameters:

  - `C`: [`J/K`] Heat capacity of element (= cp*m)
"""
@component function HeatCapacitor(; C = nothing, T = 273.15 + 20, name)
    pars = @parameters begin
        C = C, [description = "Heat capacity of element"]
    end

    systems = @named begin
        port = HeatPort()
    end

    vars = @variables begin
        T(t) = T, [guess = 273.15 + 20]
        der_T(t), [guess = 0.0]
    end

    equations = Equation[
        T ~ port.T,
        der_T ~ port.Q_flow / C,
        D(T) ~ der_T
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    ThermalConductor(; name, G)

Lumped thermal element transporting heat without storing it.

# States:

see [`Element1D`](@ref)

# Connectors:

`port_a`
`port_b`

# Parameters:

  - `G`: [`W/K`] Constant thermal conductance of material
"""
@component function ThermalConductor(; G = nothing, name)
    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d

    pars = @parameters begin
        G = G
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        Q_flow ~ G * dT
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, element1d)
end

"""
    ThermalResistor(; name, R)

Lumped thermal element transporting heat without storing it.

# States:

  - `dT`:  [`K`] Temperature difference across the component a.T - b.T
  - `Q_flow`: [`W`] Heat flow rate from port a -> port b

# Connectors:

  - `port_a`
  - `port_b`

# Parameters:

  - `R`: [`K/W`] Constant thermal resistance of material
"""
@component function ThermalResistor(; R = nothing, name)
    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d

    pars = @parameters begin
        R = R
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        dT ~ R * Q_flow
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, element1d)
end

"""
    ConvectiveConductor(; name, G)

Lumped thermal element for heat convection.

# States:

  - `dT`:  [`K`] Temperature difference across the component `solid.T` - `fluid.T`
  - `Q_flow`: [`W`] Heat flow rate from `solid` -> `fluid`

# Connectors:

  - `solid`
  - `fluid`

# Parameters:

  - `G`: [W/K] Convective thermal conductance
"""
@component function ConvectiveConductor(; G = nothing, name)
    @named convective_element1d = ConvectiveElement1D()
    @unpack Q_flow, dT = convective_element1d

    pars = @parameters begin
        G = G
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        Q_flow ~ G * dT
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, convective_element1d)
end

"""
    ConvectiveResistor(; name, R)

Lumped thermal element for heat convection.

# States:

  - `dT`: [`K`] Temperature difference across the component `solid.T` - `fluid.T`
  - `Q_flow`: [`W`] Heat flow rate from `solid` -> `fluid`

# Connectors:

  - `solid`
  - `fluid`

# Parameters:

  - `R`: [`K/W`] Constant thermal resistance of material
"""
@component function ConvectiveResistor(; R = nothing, name)
    @named convective_element1d = ConvectiveElement1D()
    @unpack Q_flow, dT = convective_element1d

    pars = @parameters begin
        R = R
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        dT ~ R * Q_flow
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, convective_element1d)
end

"""
    BodyRadiation(; name, G)

Lumped thermal element for radiation heat transfer.

# States:

  - `dT`:  [`K`] Temperature difference across the component a.T - b.T
  - `Q_flow`: [`W`] Heat flow rate from port a -> port b

# Connectors:

  - `port_a`
  - `port_b`

# Parameters:

  - `G`: [m^2] Net radiation conductance between two surfaces # Stefan-Boltzmann constant TODO: extract into physical constants module or use existing one
"""
@component function BodyRadiation(; G = nothing, name)
    sigma = 5.6703744191844294e-8 # Stefan-Boltzmann constant TODO: extract into physical constants module or use existing one

    @named element1d = Element1D()
    @unpack Q_flow, dT, port_a, port_b = element1d

    pars = @parameters begin
        G = G
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        Q_flow ~ G * sigma * (port_a.T^4 - port_b.T^4)
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, element1d)
end

"""
    ThermalCollector(; name, m = 1)

Collects `m` heat flows

This is a model to collect the heat flows from `m` heatports to one single heatport.

# States:

# Connectors:

  - `port_a1` to `port_am`
  - `port_b`

# Parameters:

  - `m`: Number of heat ports (e.g. m=2: `port_a1`, `port_a2`)
"""
@component function ThermalCollector(; m::Integer = 1, name)
    pars = @parameters begin
    end

    port_a = @named begin
        port_a[1:m] = HeatPort()
    end

    systems = @named begin
        port_b = HeatPort()
    end
    append!(systems, port_a)

    vars = @variables begin
    end

    equations = Equation[
        port_b.Q_flow + sum(k -> k.Q_flow, port_a) ~ 0,
        port_b.T ~ port_a[1].T,
        [port_a[i].T ~ port_a[i + 1].T for i in 1:(m - 1)]...
    ]

    return System(equations, t, vars, pars; name, systems)
end
