"""
    TemperatureSensor(; name)

Absolute temperature sensor in kelvin.

This is an ideal absolute temperature sensor which returns the temperature of the connected port in kelvin as an output
signal. The sensor itself has no thermal interaction with whatever it is connected to. Furthermore, no thermocouple-like
lags are associated with this sensor model.

# States:

  - `T(t)`: [`K`] Absolute temperature

# Connectors:

  - `port`
"""
@component function TemperatureSensor(; name)
    @named port = HeatPort()
    @variables T(t)
    eqs = [T ~ port.T
           port.Q_flow ~ 0]
    ODESystem(eqs, t, [T], [], systems = [port], name = name)
end

"""
    RelativeTemperatureSensor(; name)

Relative Temperature sensor.

The relative temperature `port_a.T - port_b.T` is determined between the two ports of this component and is provided as
output signal in kelvin.

# States:

  - `T(t)`: [`K`] Relative temperature `a.T - b.T`

# Connectors:

  - `port_a`
  - `port_b`
"""
@component function RelativeTemperatureSensor(; name)
    @named port_a = HeatPort()
    @named port_b = HeatPort()
    @variables T(t)
    eqs = [T ~ port_a.T - port_b.T
           port_a.Q_flow ~ 0
           port_b.Q_flow ~ 0]
    ODESystem(eqs, t, [T], [], systems = [port_a, port_b], name = name)
end

"""
    HeatFlowSensor(; name)

Heat flow rate sensor.

This model is capable of monitoring the heat flow rate flowing through this component. The sensed value of heat flow rate
is the amount that passes through this sensor while keeping the temperature drop across the sensor zero. This is an ideal
model, so it does not absorb any energy, and it has no direct effect on the thermal response of a system it is included in.
The output signal is positive, if the heat flows from `port_a` to `port_b`.

# States:

  - `Q_flow(t)`: [`W`] Heat flow from `port_a` to `port_b`

# Connectors:

  - `port_a`
  - `port_b`
"""
@component function HeatFlowSensor(; name)
    @named port_a = HeatPort()
    @named port_b = HeatPort()
    @variables Q_flow(t)
    eqs = [port_a.T ~ port_b.T
           port_a.Q_flow + port_b.Q_flow ~ 0
           Q_flow ~ port_a.Q_flow]
    ODESystem(eqs, t, [Q_flow], [], systems = [port_a, port_b], name = name)
end
