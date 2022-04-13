"""
Fixed heat flow boundary condition.

This model allows a specified amount of heat flow rate to be "injected" into a thermal system at a given port.
The constant amount of heat flow rate `Q_flow` is given as a parameter. The heat flows into the component to which
the component FixedHeatFlow is connected, if parameter `Q_flow` is positive.
"""
function FixedHeatFlow(; name, 
    Q_flow=1.0, # [W] Fixed heat flow rate at port
    T_ref=293.15, # [K] Reference temperature
    alpha=0.0, # [1/K] Temperature coefficient of heat flow rate
    )
   
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
Fixed temperature boundary condition in Kelvin.

This model defines a fixed temperature T at its port in Kelvin, i.e., it defines a fixed temperature as a boundary condition.
"""
function FixedTemperature(; name, 
    T=0.0 # [K] Fixed temperature boundary condition
    )
    @named port = HeatPort()
    pars = @parameters T=T
    eqs = [
        port.T ~ T
    ]
    ODESystem(eqs, t, [], pars; systems=[port], name=name)
end
