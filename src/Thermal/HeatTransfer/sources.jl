function FixedHeatFlow(; name, Q_flow=1.0, T₀=293.15, α=0.0)
    qflow, tem₀, alpha = Q_flow, T₀, α
    
    @parameters Q_flow T₀ α
    @named hp = HeatPort()
    
    eqs = [
        hp.Q_flow ~ -Q_flow * (1 + α*(hp.T - T₀))
    ]
    ODESystem(eqs, t, [], [Q_flow, T₀, α], systems=[hp], defaults=Dict(zip((Q_flow, T₀, α), (qflow, tem₀, alpha))), name=name)
end

function FixedTemperature(; name, T=0.0)
    tem = T

    @named hp = HeatPort()
    @parameters T

    eqs = [
        hp.T ~ T
    ]
    ODESystem(eqs, t, [], [T], systems=[hp], defaults=Dict(T => tem), name=name)
end
