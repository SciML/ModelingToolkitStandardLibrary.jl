@connector function HeatPort(; name)
    @variables T(t)=273.15 + 20.0 # [K] Temperature of the port  
    @variables Q_flow(t)=0.0 [connect = Flow] # [W] Heat flow rate at the port
    ODESystem(Equation[], t, [T, Q_flow], [], name=name)
end
Base.@doc "Port for a thermal system." HeatPort

"""
This partial model contains the basic connectors and variables to allow heat transfer models to be created that do not 
store energy. This model defines and includes equations for the temperature drop across the element, `dT`, and the heat
flow rate through the element from `port_a` to `port_b`, `Q_flow`.
"""
function Element1D(;name, 
    dT0=0.0, # [K] Temperature difference across the component a.T - b.T
    Q_flow0=0.0, # [W] Heat flow rate from port a -> port b
    )

    @named port_a = HeatPort()
    @named port_b = HeatPort()
    sts = @variables begin
        dT(t)=dT0
        Q_flow(t)=Q_flow0
    end
    eqs = [
        dT ~ port_a.T - port_b.T
        port_a.Q_flow ~ Q_flow
        port_a.Q_flow + port_b.Q_flow ~ 0
    ]
    
    return compose(ODESystem(eqs, t, sts, []; name=name), port_a, port_b)
end