using ModelingToolkitStandardLibrary.Blocks

@component function Force(; name)
    @named flange = MechanicalPort()
    @named input = RealInput()

    vars = pars = []
    eqs = [
        flange.f ~ -input.u,
    ]
    compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [flange.v => 0]),
            [flange, input])
end

@component function Position(; x_int = 0, name)
    @named flange = MechanicalPort()
    @named input = RealInput()
    pars = @parameters begin x_int = x_int end
    vars = @variables begin
        x(t) = x_int
        dx(t) = 0
    end
    eqs = [x ~ input.u
           D(x) ~ dx
           flange.v ~ dx]
    compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [flange.v => 0]),
            [flange, input])
end
