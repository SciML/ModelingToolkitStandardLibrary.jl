"""
    PressureSource(; name, dp)

Provides a pressure difference between port `a` and port `b`.

# Connectors:
- `a` left port [`HydraulicPort`](@ref)
- `b` right port [`HydraulicPort`](@ref)

# Parameters:
- `dp`: [`Pa`] Pressure difference `a.p ~ b.p + dp`
"""
function PressureSource(; name, dp)
    @named a = HydraulicPort()
    @named b = HydraulicPort()
    @parameters dp = dp
    eqs = [
        a.p ~ b.p + dp
        a.m_flow + b.m_flow ~ 0
    ]
    return ODESystem(eqs, t, [], [dp]; systems = [a, b], name = name)
end

"""
    MassFlowRateSource(; name, m_flow)

Provides a mass flow. A positive value causes liquid to flow from `a` to `b`.

# Connectors:
- `a` left port [`HydraulicPort`](@ref)
- `b` right port [`HydraulicPort`](@ref)

# Parameters:
- `m_flow`: [`kg/s`] Prescribed mass flow rate
"""
function MassFlowRateSource(; name, m_flow)
    @named a = HydraulicPort()
    @named b = HydraulicPort()
    @parameters m_flow = m_flow
    eqs = [
        a.m_flow + b.m_flow ~ 0
        a.m_flow ~ m_flow
    ]
    return ODESystem(eqs, t, [], [m_flow]; systems = [a, b], name = name)
end
