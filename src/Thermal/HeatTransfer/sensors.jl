"""
Absolute temperature sensor in Kelvin.
"""
function TemperatureSensor(; name)
    @named a = HeatPort()
    @variables T(t) # [K] Absolute temperature
    
    eqs = [
        T ~ a.T
        a.Q_flow ~ 0
    ]
    ODESystem(eqs, t, [T], [], systems=[a], name=name)
end

"""
Relative Temperature sensor.
"""
function RelativeTemperatureSensor(; name)
    @named a = HeatPort()
    @named b = HeatPort()
    @variables T(t) # [K] Relative temperature a.T - b.T
    
    eqs = [
        T ~ a.T - b.T
        a.Q_flow ~ 0
        b.Q_flow ~ 0
    ]
    ODESystem(eqs, t, [T], [], systems=[a, b], name=name)
end

"""
Heat flow rate sensor.
"""
function HeatFlowSensor(; name)
    @named a = HeatPort()
    @named b = HeatPort()
    @variables Q_flow(t) # [W] Heat flow from port a to port b
    
    eqs = [
        a.T ~ b.T
        a.Q_flow + b.Q_flow ~ 0
        Q_flow ~ a.Q_flow
    ]
    ODESystem(eqs, t, [Q_flow], [], systems=[a, b], name=name)
end
