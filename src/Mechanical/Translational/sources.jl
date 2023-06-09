"""
    Force(; name)

Linear 1D force input source

# Connectors:

  - `flange`: 1-dim. translational flange
  - `f`: real input 
"""
@component function Force(; name)
    systems = @named begin
        flange = MechanicalPort()
        f = RealInput()
    end

    vars = pars = []
    eqs = [
        flange.f ~ -f.u,
    ]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [flange.v => 0])
end

"""
    Position(; s_0 = 0, name)

Linear 1D position input source

# Parameters:

- `s_0`: [m] initial value of absolute position

# Connectors:

  - `flange`: 1-dim. translational flange
  - `s`: real input 
"""
@component function Position(; s_0 = 0, name)
    systems = @named begin
        flange = MechanicalPort()
        s = RealInput()
    end

    pars = @parameters s_0 = s_0
    vars = @variables x(t) = s_0

    eqs = [D(x) ~ flange.v
        s.u ~ x]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [flange.v => 0, s.u => s_0])
end
