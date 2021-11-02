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

function PotentialSensor(; name)
    @named p = Pin()
    @variables phi(t)
    eqs = [
        p.i ~ 0
        phi ~ p.v
    ]
    ODESystem(eqs, t, [phi], [], systems=[p], defaults=Dict(phi => 1.0), name=name)
end

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
