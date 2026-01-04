"""
    ForceSensor(; name)

Linear 1D force sensor, measures the force between two flanges.

# Connectors:

- `flange`: 1-dim. translational flange
- `output`: real output
"""
@component function ForceSensor(; name)
    pars = @parameters begin
    end

    systems = @named begin
        flange_a = MechanicalPort()
        flange_b = MechanicalPort()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        flange_a.v ~ flange_b.v,
        flange_a.f + flange_b.f ~ 0.0,
        output.u ~ flange_a.f,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    PositionSensor(; s = 0, name)

Linear 1D position sensor.

# States:

- `s`: [m] absolute position (with initial value of 0.0)

# Connectors:

- `flange`: 1-dim. translational flange
- `output`: real output
"""
@component function PositionSensor(; name, s = nothing)
    pars = @parameters begin
    end

    systems = @named begin
        flange = MechanicalPort()
        output = RealOutput()
    end

    vars = @variables begin
        s(t) = s
    end

    equations = Equation[
        D(s) ~ flange.v,
        output.u ~ s,
        flange.f ~ 0.0,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    AccelerationSensor(; name)

Linear 1D acceleration sensor.

# States:

- `a`: [m/s^2] measured acceleration

# Connectors:

- `flange`: 1-dim. translational flange
- `output`: real output
"""
@component function AccelerationSensor(; name, a = nothing)
    pars = @parameters begin
    end

    systems = @named begin
        flange = MechanicalPort()
        output = RealOutput()
    end

    vars = @variables begin
        a(t) = a
    end

    equations = Equation[
        a ~ D(flange.v),
        output.u ~ a,
        flange.f ~ 0.0,
    ]

    return System(equations, t, vars, pars; name, systems)
end
