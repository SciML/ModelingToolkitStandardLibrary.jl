"""
    FixedHeatFlow(; name, Q_flow=1.0, T_ref=293.15, alpha=0.0)

Fixed heat flow boundary condition.

This model allows a specified amount of heat flow rate to be "injected" into a thermal system at a given port.
The constant amount of heat flow rate `Q_flow` is given as a parameter. The heat flows into the component to which
the component FixedHeatFlow is connected, if parameter `Q_flow` is positive.

# Connectors:
- `port`

# Parameters:
- `Q_flow`: [W] Fixed heat flow rate at port
- `T_ref`: [K] Reference temperature
- `alpha`: [1/K] Temperature coefficient of heat flow rate
"""
function FixedHeatFlow(; name, Q_flow=1.0, T_ref=293.15, alpha=0.0)
    pars = @parameters begin 
        Q_flow=Q_flow 
        T_ref=T_ref 
        alpha=alpha
    end
    @named port = HeatPort()
    
    eqs = [
        port.Q_flow ~ -Q_flow * (1 + alpha * (port.T - T_ref))
    ]
    ODESystem(eqs, t, [], pars; systems=[port], name=name)
end

"""
    FixedTemperature(; name, T)

Fixed temperature boundary condition in kelvin.

This model defines a fixed temperature `T` at its port in kelvin, i.e., it defines a fixed temperature as a boundary condition.

# Connectors:
- `port`

# Parameters:
- `T`: [K] Fixed temperature boundary condition
"""
function FixedTemperature(; name, T)
    @named port = HeatPort()
    pars = @parameters T=T
    eqs = [
        port.T ~ T
    ]
    ODESystem(eqs, t, [], pars; systems=[port], name=name)
end


"""
    PrescribedHeatFlow(; name, Q_flow=1.0, T_ref=293.15, alpha=0.0)

Prescribed heat flow boundary condition.

This model allows a specified amount of heat flow rate to be "injected" into a thermal system at a given port. 
The amount of heat is given by the input signal `Q_flow` into the model. The heat flows into the component to which 
the component `PrescribedHeatFlow` is connected, if the input signal is positive.
If parameter alpha is > 0, the heat flow is multiplied by `1 + alpha*(port.T - T_ref`) in order to simulate temperature 
dependent losses (which are given an reference temperature T_ref).

# Connectors:
- `port`
- `Q_flow` Input for the heat flow

# Parameters:
- `T_ref`: [K] Reference temperature
- `alpha`: [1/K] Temperature coefficient of heat flow rate
"""
function PrescribedHeatFlow(; name, T_ref=293.15, alpha=0.0)
    pars = @parameters begin 
        Q_flow=Q_flow 
        T_ref=T_ref 
        alpha=alpha
    end
    @named port = HeatPort()
    @named Q_flow = RealInput()
    
    eqs = [
        port.Q_flow ~ -Q_flow.u * (1 + alpha * (port.T - T_ref))
    ]
    ODESystem(eqs, t, [], pars; systems=[port, Q_flow], name=name)
end

"""
    PrescribedTemperature(; name, T)

This model represents a variable temperature boundary condition. 

The temperature in kelvin is given as input signal `T` to the model. The effect is that an instance of 
this model acts as an infinite reservoir able to absorb or generate as much energy as required to keep 
the temperature at the specified value.

# Connectors:
- `port`
- `T` input for the temperature
"""
function PrescribedTemperature(; name)
    @named port = HeatPort()
    @named T = RealInput()
    eqs = [
        port.T ~ T.u
    ]
    ODESystem(eqs, t, [], pars; systems=[port, T], name=name)
end