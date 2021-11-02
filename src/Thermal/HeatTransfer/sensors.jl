function TemperatureSensor(; name)
    @named hp = HeatPort()
    @variables T(t)
    
    eqs = [
        T ~ hp.T
        hp.Q_flow ~ 0
    ]
    ODESystem(eqs, t, [T], [], systems=[hp], name=name)
end

function RelativeTemperatureSensor(; name)
    @named hp1 = HeatPort()
    @named hp2 = HeatPort()
    @variables T(t)
    
    eqs = [
        T ~ hp1.T - hp2.T
        hp1.Q_flow ~ 0
        hp2.Q_flow ~ 0
    ]
    ODESystem(eqs, t, [T], [], systems=[hp1, hp2], name=name)
end

function HeatFlowSensor(; name)
    @named hp1 = HeatPort()
    @named hp2 = HeatPort()
    @variables Q_flow(t)
    
    eqs = [
        hp1.T ~ hp2.T
        hp1.Q_flow + hp2.Q_flow ~ 0
        Q_flow ~ hp1.Q_flow
    ]
    ODESystem(eqs, t, [Q_flow], [], systems=[hp1, hp2], name=name)
end
