function PositionSensor(; name, s0 = 0.0)
    @named port = Port()

    pars = @parameters s0 = s0
    vars = @variables s(t) = s0
    eqs = [D(s) ~ port.v
           0 ~ p.f]
    compose(ODESystem(eqs, t, vars, pars; name = name), port)
end
