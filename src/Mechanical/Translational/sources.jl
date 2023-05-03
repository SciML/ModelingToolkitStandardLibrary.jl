"""
    Force(; name)

Linear 1D force input source

# Connectors:

  - `flange`: 1-dim. translational flange
  - `input`: real input 
"""
@component function Force(; name)
    systems = @named begin
        flange = MechanicalPort()
        input = RealInput()
    end

    vars = pars = []
    eqs = [
        flange.f ~ -input.u,
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
  - `input`: real input 
"""
@component function Position(; s_0 = 0, name)
    systems = @named begin
        flange = MechanicalPort()
        input = RealInput()
    end

    pars = @parameters s_0 = s_0
    vars = @variables s(t) = s_0

    eqs = [D(s) ~ flange.v
           input.u ~ s]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [flange.v => 0, input.u => s_0])
end
