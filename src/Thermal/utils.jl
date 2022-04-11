@connector function HeatPort(; name)
    @variables T(t)=273.15, Q_flow(t)=0.0 [connect=Flow] # Temperature and Heat-flow-rate
    ODESystem(Equation[], t, [T, Q_flow], [], name=name)
end