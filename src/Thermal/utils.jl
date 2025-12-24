@connector function HeatPort(; name, T_guess = 273.15 + 20, Q_flow_guess = 0.0, T = nothing, Q_flow = nothing)
    pars = @parameters begin
        T_guess = 273.15 + 20
        Q_flow_guess = 0.0
    end

    vars = @variables begin
        T(t) = T, [guess = T_guess]
        Q_flow(t) = Q_flow, [guess = Q_flow_guess, connect = Flow]
    end
    System(Equation[], t, vars, pars; name)
end
Base.@doc """
    HeatPort(; T = nothing, T_guess = 273.15 + 20, Q_flow = nothing, Q_flow_guess = 0.0, name)

Port for a thermal system.
# Parameters: 
- `T_guess`: [K] Initial guess for the temperature of the port (set to 273.15 + 20).
- `Q_flow_guess`: [W] Initial guess for the heat flow rate at the port (set to 0.0).

# States:
- `T`: [K] Temperature of the port. Guess set to `T_guess`. Passing a value for `T` will set its default.
- `Q_flow`: [W] Heat flow rate at the port. Guess set to `Q_flow_guess`. Passing a value for `Q_flow` will set its default.
""" HeatPort

"""
    Element1D(; name, dT = 0.0, Q_flow = 0.0)

This partial model contains the basic connectors and variables to allow heat transfer models to be created that do not
store energy. This model defines and includes equations for the temperature drop across the element, `dT`, and the heat
flow rate through the element from `port_a` to `port_b`, `Q_flow`.

# States:

  - `dT`:  [`K`] Temperature difference across the component a.T - b.T. It accepts an initial value, which defaults to 0.0.
  - `Q_flow`: [`W`] Heat flow rate from port a -> port b. It accepts an initial value, which defaults to 0.0.

# Connectors:

`port_a`
`port_b`
"""
@component function Element1D(; dT = 0.0, Q_flow = 0.0, name)
    pars = @parameters begin
    end

    systems = @named begin
        port_a = HeatPort()
        port_b = HeatPort()
    end

    vars = @variables begin
        dT(t) = dT, [guess = 0.0]
        Q_flow(t) = Q_flow, [guess = 0.0]
    end

    equations = Equation[
        dT ~ port_a.T - port_b.T,
        port_a.Q_flow ~ Q_flow,
        port_a.Q_flow + port_b.Q_flow ~ 0
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    ConvectiveElement1D(; name, dT = 0.0, Q_flow = 0.0)

This partial model contains the basic connectors and variables to allow heat
transfer models to be created that do not store energy. This model defines and
includes equations for the temperature drop across the element, `dT`, and the heat
flow rate through the element from `solid` to `fluid`, `Q_flow`.

# States:

  - `dT`:  [`K`] Temperature difference across the component `solid.T` - `fluid.T`. It accepts an initial value, which defaults to 0.0.
  - `Q_flow`: [`W`] Heat flow rate from `solid` -> `fluid`. It accepts an initial value, which defaults to 0.0.

# Connectors:

`solid`
`fluid`
"""
@component function ConvectiveElement1D(; dT = 0.0, Q_flow = 0.0, name)
    pars = @parameters begin
    end

    systems = @named begin
        solid = HeatPort()
        fluid = HeatPort()
    end

    vars = @variables begin
        dT(t) = dT, [guess = 0.0]
        Q_flow(t) = Q_flow, [guess = 0.0]
    end

    equations = Equation[
        dT ~ solid.T - fluid.T,
        solid.Q_flow ~ Q_flow,
        solid.Q_flow + fluid.Q_flow ~ 0
    ]

    return System(equations, t, vars, pars; name, systems)
end
