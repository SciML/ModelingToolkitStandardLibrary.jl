"""
```julia
function CurrentSensor(; name)
```

Creates a circuit component that measures the current flowing through it. Analogous to
an ideal ammeter.

# States
- `i(t)`: [`A`]
  Current through the sensor

# Connectors
- `p`
 Positive pin
- `n`
  Negative pin
"""
function CurrentSensor(; name)
    @named p = Pin()
    @named n = Pin()
    @variables i(t)=1.0
    eqs = [
        p.v ~ n.v
        i ~ p.i
        i ~ -n.i
    ]
    ODESystem(eqs, t, [i], [], systems=[p, n]; name=name)
end

"""
```julia
function PotentialSensor(; name)
```

Creates a circuit component which measures the potential at a pin.

# States
- `phi(t)`: [`V`]
  The potential at this point

# Connectors
- `p`
  Pin at which potential is to be measured
"""
function PotentialSensor(; name)
    @named p = Pin()
    @variables phi(t)=1.0
    eqs = [
        p.i ~ 0
        phi ~ p.v
    ]
    ODESystem(eqs, t, [phi], [], systems=[p]; name=name)
end

"""
```julia
function VoltageSensor(; name)
```

Creates a circuit component that measures the voltage across it. Analogous to
an ideal voltmeter.

# States
- `v(t)`: [`V`]
  The voltage across this component

# Connectors
- `p`
  Positive pin
- `n`
  Negative pin
"""
function VoltageSensor(; name)
    @named p = Pin()
    @named n = Pin()
    @variables v(t)=1.0
    eqs = [
        p.i ~ 0
        n.i ~ 0
        v ~ p.v - n.v
    ]
    ODESystem(eqs, t, [v], []; systems=[p, n], name=name)
end

"""
```julia
function PowerSensor(; name)
```

Combines a [`VoltageSensor`](@ref) and a [`CurrentSensor`](@ref) to measure the power being
consumed by a circuit.

# States
- `power(t)`: [`W`]
  The power being consumed, given by the product of voltage and current.

# Connectors
- `pc`
  Corresponds to the `p` pin of the [`CurrentSensor`](@ref)
- `nc`
  Corresponds to the `n` pin of the [`CurrentSensor`](@ref)
- `pv`
  Corresponds to the `p` pin of the [`VoltageSensor`](@ref)
- `nv`
  Corresponds to the `n` pin of the [`VoltageSensor`](@ref)
"""
function PowerSensor(; name)
    @named pc = Pin()
    @named nc = Pin()
    @named pv = Pin()
    @named nv = Pin()
    @named voltage_sensor = VoltageSensor()
    @named current_sensor = CurrentSensor()
    @variables power(t)=1.0
    eqs = [
        connect(voltage_sensor.p, pv)
        connect(voltage_sensor.n, nv)
        connect(current_sensor.p, pc)
        connect(current_sensor.n, nc)  
        power ~ current_sensor.i * voltage_sensor.v
    ]
    ODESystem(eqs, t, [power], []; systems=[pc, nc, pv, nv, voltage_sensor, current_sensor], name=name)
end

"""
```julia
function MultiSensor(; name)
```

Combines a [`VoltageSensor`](@ref) and a [`CurrentSensor`](@ref).

# States
- `v(t)`: [`V`]
  The voltage across the [`VoltageSensor`](@ref)
- `i(t)`: [`A`]
  The current across the [`CurrentSensor`](@ref)

# Connectors
- `pc`
  Corresponds to the `p` pin of the [`CurrentSensor`](@ref)
- `nc`
  Corresponds to the `n` pin of the [`CurrentSensor`](@ref)
- `pv`
  Corresponds to the `p` pin of the [`VoltageSensor`](@ref)
- `nv`
  Corresponds to the `n` pin of the [`VoltageSensor`](@ref)
"""
function MultiSensor(; name)
    @named pc = Pin()
    @named nc = Pin()
    @named pv = Pin()
    @named nv = Pin()
    @named voltage_sensor = VoltageSensor()
    @named current_sensor = CurrentSensor()
    sts = @variables begin
        i(t)=1.0
        v(t)=1.0
    end
    eqs = [
        connect(voltage_sensor.p, pv)
        connect(voltage_sensor.n, nv)
        connect(current_sensor.p, pc)
        connect(current_sensor.n, nc)
        i ~ current_sensor.i
        v ~ voltage_sensor.v
    ]
    ODESystem(eqs, t, sts, []; systems=[pc, nc, pv, nv, voltage_sensor, current_sensor], name=name)
end
