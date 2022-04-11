"""
Fixed heat flow boundary condition.
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
    @named b = HeatPort()
    
    eqs = [
        b.Q_flow ~ -Q_flow * (1 + alpha * (b.T - T_ref))
    ]
    ODESystem(eqs, t, [], pars; systems=[b], name=name)
end

"""
Fixed temperature boundary condition in Kelvin.
"""
function FixedTemperature(; name, 
    T=0.0 # [K] Fixed temperature boundary condition
    )
    @named b = HeatPort()
    pars = @parameters T=T
    eqs = [
        b.T ~ T
    ]
    ODESystem(eqs, t, [], pars; systems=[b], name=name)
end
