@connector function HeatPort(; name, T_start = 273.15 + 20.0, Q_flow_start = 0.0)
    @variables T(t) = T_start
    @variables Q_flow(t)=Q_flow_start [connect = Flow]
    ODESystem(Equation[], t, [T, Q_flow], [], name = name)
end
Base.@doc """
    HeatPort(; name, T_start=273.15 + 20.0, Q_flow_start=0.0)

Port for a thermal system.

# States:
- `T`: [K] Temperature of the port
- `Q_flow`: [W] Heat flow rate at the port

# Parameters:
- `T_start`: [K] Temperature of the port
- `Q_flow_start`: [W] Heat flow rate at the port
""" HeatPort

"""
    Element1D(; name, dT_start = 0.0, Q_flow_start = 0.0)

This partial model contains the basic connectors and variables to allow heat transfer models to be created that do not
store energy. This model defines and includes equations for the temperature drop across the element, `dT`, and the heat
flow rate through the element from `port_a` to `port_b`, `Q_flow`.

# States:

  - `dT`:  [`K`] Temperature difference across the component a.T - b.T
  - `Q_flow`: [`W`] Heat flow rate from port a -> port b

# Connectors:

`port_a`
`port_b`

# Parameters:

  - `dT_start`:  [K] Initial temperature difference across the component a.T - b.T
  - `Q_flow_start`: [W] Initial heat flow rate from port a -> port b
"""
@component function Element1D(; name, dT_start = 0.0, Q_flow_start = 0.0)
    @named port_a = HeatPort()
    @named port_b = HeatPort()
    sts = @variables begin
        dT(t) = dT_start
        Q_flow(t) = Q_flow_start
    end
    eqs = [dT ~ port_a.T - port_b.T
        port_a.Q_flow ~ Q_flow
        port_a.Q_flow + port_b.Q_flow ~ 0]

    return compose(ODESystem(eqs, t, sts, []; name = name), port_a, port_b)
end

"""
    ConvectiveElement1D(; name, dT_start = 0.0, Q_flow_start = 0.0)

This partial model contains the basic connectors and variables to allow heat
transfer models to be created that do not store energy. This model defines and
includes equations for the temperature drop across the element, `dT`, and the heat
flow rate through the element from `solid` to `fluid`, `Q_flow`.

# States:

  - `dT`:  [`K`] Temperature difference across the component `solid.T` - `fluid.T`
  - `Q_flow`: [`W`] Heat flow rate from `solid` -> `fluid`

# Connectors:

`solid`
`fluid`

# Parameters:

  - `dT_start`:  [K] Initial temperature difference across the component `solid.T` - `fluid.T`
  - `Q_flow_start`: [W] Initial heat flow rate from `solid` -> `fluid`
"""
@component function ConvectiveElement1D(; name, dT_start = 0.0, Q_flow_start = 0.0)
    @named solid = HeatPort()
    @named fluid = HeatPort()
    sts = @variables begin
        dT(t) = dT_start
        Q_flow(t) = Q_flow_start
    end
    eqs = [dT ~ solid.T - fluid.T
        solid.Q_flow ~ Q_flow
        solid.Q_flow + fluid.Q_flow ~ 0]

    return compose(ODESystem(eqs, t, sts, []; name = name), solid, fluid)
end
