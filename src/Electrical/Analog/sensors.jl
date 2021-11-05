"""
```julia
function CurrentSensor(; name)
```

Creates a circuit component that measures the current flowing through it. Analogous to
an ideal ammeter.

# States
- `i(t)`
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
    @variables i(t)
    eqs = [
        p.v ~ n.v
        i ~ p.i
        i ~ -n.i
    ]
    ODESystem(eqs, t, [i], [], systems=[p, n], defaults=Dict(i => 1.0), name=name)
end

"""
```julia
function PotentialSensor(; name)
```

Creates a circuit component which measures the potential at a pin.

# States
- `phi(t)`
  The potential at this point

# Connectors
- `p`
  Pin at which potential is to be measured
"""
function PotentialSensor(; name)
    @named p = Pin()
    @variables phi(t)
    eqs = [
        p.i ~ 0
        phi ~ p.v
    ]
    ODESystem(eqs, t, [phi], [], systems=[p], defaults=Dict(phi => 1.0), name=name)
end

"""
```julia
function VoltageSensor(; name)
```

Creates a circuit component that measures the voltage across it. Analogous to
an ideal voltmeter.

# States
- `v(t)`
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
    @variables v(t)
    eqs = [
        p.i ~ 0
        n.i ~ 0
        v ~ p.v - n.v
    ]
    ODESystem(eqs, t, [v], [], systems=[p, n], defaults=Dict(v => 1.0), name=name)
end

"""
```julia
function PowerSensor(; name)
```

Combines a [`VoltageSensor`](@ref) and a [`CurrentSensor`](@ref) to measure the power being
consumed by a circuit.

# States
- `power(t)`
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
    @variables power(t)
    eqs = [
        connect(voltage_sensor.p, pv)
        connect(voltage_sensor.n, nv)
        connect(current_sensor.p, pc)
        connect(current_sensor.n, nc)  
        power ~ current_sensor.i * voltage_sensor.v
    ]
    ODESystem(eqs, t, [power], [], systems=[pc, nc, pv, nv, voltage_sensor, current_sensor], defaults=Dict(power => 1.0), name=name)
end

"""
```julia
function MultiSensor(; name)
```

Combines a [`VoltageSensor`](@ref) and a [`CurrentSensor`](@ref).

# States
- `v(t)`
  The voltage across the [`VoltageSensor`](@ref)
- `i(t)`
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
    @variables i(t) v(t)
    eqs = [
        connect(voltage_sensor.p, pv)
        connect(voltage_sensor.n, nv)
        connect(current_sensor.p, pc)
        connect(current_sensor.n, nc)
        i ~ current_sensor.i
        v ~ voltage_sensor.v
    ]
    ODESystem(eqs, t, [i, v], [], systems=[pc, nc, pv, nv, voltage_sensor, current_sensor], defaults=Dict(i => 1.0, v => 1.0), name=name)
end
