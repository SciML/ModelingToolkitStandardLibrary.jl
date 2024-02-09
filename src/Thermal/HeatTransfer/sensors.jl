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
@mtkmodel TemperatureSensor begin
    @components begin
        port = HeatPort()
    end
    @variables begin
        T(t), [description = "Absolute temperature", unit = u"K"]
    end
    @equations begin
        T ~ port.T
        port.Q_flow ~ 0
    end
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
@mtkmodel RelativeTemperatureSensor begin
    @components begin
        port_a = HeatPort()
        port_b = HeatPort()
    end
    @variables begin
        T(t), [description = "Relative temperature", unit = u"K"]
    end
    @equations begin
        T ~ port_a.T - port_b.T
        port_a.Q_flow ~ 0
        port_b.Q_flow ~ 0
    end
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
@mtkmodel HeatFlowSensor begin
    @components begin
        port_a = HeatPort()
        port_b = HeatPort()
    end
    @variables begin
        Q_flow(t), [connect = Flow, description = "Heat flow rate", unit = u"W"]
    end
    @equations begin
        port_a.T ~ port_b.T
        port_a.Q_flow + port_b.Q_flow ~ 0
        Q_flow ~ port_a.Q_flow
    end
end
