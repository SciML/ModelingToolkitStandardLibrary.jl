"""
    FixedHeatFlow(; name, Q_flow=1.0, T_ref=293.15, alpha=0.0)

Fixed heat flow boundary condition.

This model allows a specified amount of heat flow rate to be "injected" into a thermal system at a given port.
The constant amount of heat flow rate `Q_flow` is given as a parameter. The heat flows into the component to which
the component FixedHeatFlow is connected, if parameter `Q_flow` is positive.

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

This model defines a fixed temperature T at its port in kelvin, i.e., it defines a fixed temperature as a boundary condition.

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
