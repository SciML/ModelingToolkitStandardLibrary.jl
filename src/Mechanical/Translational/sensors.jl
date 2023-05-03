"""
    ForceSensor(; name)

Linear 1D force input sensor.

# Connectors:

  - `flange`: 1-dim. translational flange
  - `output`: real output
"""
@component function ForceSensor(; name)
    systems = @named begin
        flange = MechanicalPort()
        output = RealOutput()
    end

    vars = pars = []
    eqs = [
        flange.f ~ -output.u,
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    PositionSensor(; s_0 = 0, name)

Linear 1D position input sensor.

# Parameters:

- `s_0`: [m] initial value of absolute position

# Connectors:

  - `flange`: 1-dim. translational flange
  - `output`: real output
"""
@component function PositionSensor(; s_0 = 0, name)
    systems = @named begin
        flange = MechanicalPort()
        output = RealOutput()
    end

    pars = @parameters s_0 = s_0
    vars = @variables s(t) = s_0

    eqs = [D(s) ~ flange.v
           output.u ~ s
           flange.f ~ 0]

    ODESystem(eqs, t, vars, pars; name, systems,
              defaults = [flange.v => 0, output.u => s_0])
end
