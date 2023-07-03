"""
    ForceSensor(; name)

Linear 1D force input sensor.

# Connectors:

  - `flange`: 1-dim. translational flange
  - `output`: real output
"""
@mtkmodel ForceSensor begin
    @components begin
        flange = MechanicalPort()
        output = RealOutput()
    end

    @equations begin
        flange.f ~ -output.u
    end
end

"""
    PositionSensor(; s = 0, name)

Linear 1D position input sensor.

# States:

- `s`: [m] absolute position (with initial value of 0.0)

# Connectors:

  - `flange`: 1-dim. translational flange
  - `output`: real output
"""
@mtkmodel PositionSensor begin
    @components begin
        flange = MechanicalPort()
        output = RealOutput()
    end

    @variables begin
        s(t) = 0.0
    end

    @equations begin
        D(s) ~ flange.v
        output.u ~ s
        flange.f ~ 0.0
    end

end
