@connector function HeatPort(; name)
    @variables T(t)=273.15 + 20.0 # [K] Temperature of the port  
    @variables Q_flow(t)=0.0 [connect = Flow] # [W] Heat flow rate at the port
    ODESystem(Equation[], t, [T, Q_flow], [], name=name)
end
Base.@doc "Port for a thermal system." HeatPort

function Element1D(;name, 
    dT0=0.0, # [K] Temperature difference across the component a.T - b.T
    Q_flow0=0.0, # [W] Heat flow rate from port a -> port b
    )

    @named a = HeatPort()
    @named b = HeatPort()
    sts = @variables begin
        dT(t)=dT0
        Q_flow(t)=Q_flow0
    end
    eqs = [
        dT ~ a.T - b.T
        a.Q_flow ~ Q_flow
        a.Q_flow + b.Q_flow ~ 0
    ]
    
    return compose(ODESystem(eqs, t, sts, []; name=name), a, b)
end