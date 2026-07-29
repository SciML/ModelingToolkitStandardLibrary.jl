"""
    HeatPort(; name, T = nothing, T_guess = 293.15, Q_flow = nothing, Q_flow_guess = 0.0)

Port for a thermal system.

# Keyword Arguments
- `T_guess`: [`K`] Initial temperature guess.
- `Q_flow_guess`: [`W`] Initial heat-flow guess.
- `T`: Default value for the temperature state.
- `Q_flow`: Default value for the heat-flow state.
"""
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

"""
    Element1D(; name, dT_guess = 0.0, Q_flow_guess = 0.0)

This partial model contains the basic connectors and variables to allow heat transfer models to be created that do not
store energy. This model defines and includes equations for the temperature drop across the element, `dT`, and the heat
flow rate through the element from `port_a` to `port_b`, `Q_flow`.

# States:

  - `dT`:  [`K`] Temperature difference across the component a.T - b.T (algebraically constrained).
  - `Q_flow`: [`W`] Heat flow rate from port a -> port b (algebraically constrained).

# Connectors:

`port_a`
`port_b`
"""
@component function Element1D(; name, dT_guess = 0.0, Q_flow_guess = 0.0)
    pars = @parameters begin
    end

    systems = @named begin
        port_a = HeatPort()
        port_b = HeatPort()
    end

    vars = @variables begin
        dT(t), [guess = dT_guess]
        Q_flow(t), [guess = Q_flow_guess]
    end

    equations = Equation[
        dT ~ port_a.T - port_b.T,
        port_a.Q_flow ~ Q_flow,
        port_a.Q_flow + port_b.Q_flow ~ 0,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    ConvectiveElement1D(; name, dT_guess = 0.0, Q_flow_guess = 0.0)

This partial model contains the basic connectors and variables to allow heat
transfer models to be created that do not store energy. This model defines and
includes equations for the temperature drop across the element, `dT`, and the heat
flow rate through the element from `solid` to `fluid`, `Q_flow`.

# States:

  - `dT`:  [`K`] Temperature difference across the component `solid.T` - `fluid.T` (algebraically constrained).
  - `Q_flow`: [`W`] Heat flow rate from `solid` -> `fluid` (algebraically constrained).

# Connectors:

`solid`
`fluid`
"""
@component function ConvectiveElement1D(; name, dT_guess = 0.0, Q_flow_guess = 0.0)
    pars = @parameters begin
    end

    systems = @named begin
        solid = HeatPort()
        fluid = HeatPort()
    end

    vars = @variables begin
        dT(t), [guess = dT_guess]
        Q_flow(t), [guess = Q_flow_guess]
    end

    equations = Equation[
        dT ~ solid.T - fluid.T,
        solid.Q_flow ~ Q_flow,
        solid.Q_flow + fluid.Q_flow ~ 0,
    ]

    return System(equations, t, vars, pars; name, systems)
end
