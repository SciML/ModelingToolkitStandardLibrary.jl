"""
    MassFlow(; name, p_int)

Hydraulic mass flow input source

# Connectors:

  - `port`: hydraulic port
  - `dm`: real input 
"""
@component function MassFlow(; name, p_int)
    pars = @parameters p_int = p_int

    systems = @named begin
        port = HydraulicPort(; p_int)
        dm = RealInput()
    end

    vars = []
    eqs = [
        port.dm ~ -dm.u
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    FixedPressure(; p, name)

Fixed pressure source

# Parameters:
- `p`: [Pa] set pressure (set by `p` argument)

# Connectors:
- `port`: hydraulic port
"""
@component function FixedPressure(; p, name)
    pars = @parameters begin
        p = p
    end

    vars = []

    systems = @named begin
        port = HydraulicPort()
    end

    eqs = [
        port.p ~ p
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end
@deprecate Source FixedPressure

"""
    Pressure(; p_int, name)

input pressure source

# Parameters:
- `p_int`: [Pa] initial pressure (set by `p_int` argument)

# Connectors:
- `port`: hydraulic port
- `p`: real input 
"""
@component function Pressure(; name)
    pars = []
    vars = []

    systems = @named begin
        port = HydraulicPort()
        p = RealInput()
    end

    eqs = [
        port.p ~ p.u
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end
@deprecate InputSource Pressure
