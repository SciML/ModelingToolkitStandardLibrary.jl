"""
    MassFlow(; name)

Hydraulic mass flow input source

# Connectors:

  - `port`: hydraulic port
  - `dm`: real input
"""
@component function MassFlow(; name)
    pars = @parameters begin
    end

    @named port = HydraulicPort()
    @named dm = RealInput()
    systems = [port, dm]

    vars = @variables begin
    end

    equations = Equation[
        port.dm ~ -dm.u,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    FixedPressure(; p, name)

Fixed pressure source

# Parameters:
- `p`: [Pa] set pressure (set by `p` argument)

# Connectors:
- `port`: hydraulic port
"""
@component function FixedPressure(; name, p = nothing)
    pars = @parameters begin
        p = p
    end

    @named port = HydraulicPort()
    systems = [port]

    vars = @variables begin
    end

    equations = Equation[
        port.p ~ p,
    ]

    return System(equations, t, vars, pars; name, systems)
end
@deprecate Source FixedPressure

"""
    Pressure(; name)

input pressure source

# Connectors:
- `port`: hydraulic port
- `p`: real input
"""
@component function Pressure(; name)
    pars = @parameters begin
    end

    @named port = HydraulicPort()
    @named p = RealInput()
    systems = [port, p]

    vars = @variables begin
    end

    equations = Equation[
        port.p ~ p.u,
    ]

    return System(equations, t, vars, pars; name, systems)
end
@deprecate InputSource Pressure
