@component function PositionSensor(; name, s_0 = 0.0)
    @named flange = MechanicalPort()

    pars = @parameters s_0 = s_0
    vars = @variables s(t) = s_0
    eqs = [D(s) ~ flange.v
           0 ~ p.f]
    compose(ODESystem(eqs, t, vars, pars; name = name), flange)
end
