using ModelingToolkitStandardLibrary.Blocks

function Force(; name)
    @named flange = MechanicalPort()
    @named f = RealInput()

    vars = pars = []
    eqs = [
        flange.f ~ -f.u,
    ]
    compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [flange.v => 0]),
            [flange, f])
end
