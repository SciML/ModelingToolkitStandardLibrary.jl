@connector function HeatPort(; name)
    sts = @variables begin 
        T(t)                        # Temperature in [K]
        Q_flow(t), [connect=Flow]   # Heat flow rate in [W]
    end
    ODESystem(Equation[], t, sts, [], name=name)
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
        b.Q_flow ~ -Q_flow
    ]
    
    return compose(ODESystem(eqs, t, sts, []; name=name), a, b)
end