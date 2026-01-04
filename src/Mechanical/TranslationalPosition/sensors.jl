"""
    ForceSensor(; name)

Linear 1D force sensor, measures the force between two flanges.

# Connectors:

- `flange_a`: 1-dim. translational flange
- `flange_b`: 1-dim. translational flange
- `output`: real output
"""
@component function ForceSensor(; name)
    pars = @parameters begin
    end

    systems = @named begin
        flange_a = Flange()
        flange_b = Flange()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        flange_a.s ~ flange_b.s,
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
@component function PositionSensor(; name)
    pars = @parameters begin
    end

    systems = @named begin
        flange = Flange()
        output = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u ~ flange.s,
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
        flange = Flange()
        output = RealOutput()
    end

    vars = @variables begin
        a(t) = a
    end

    equations = Equation[
        a ~ D(D(flange.s)),
        output.u ~ a,
        flange.f ~ 0.0,
    ]

    return System(equations, t, vars, pars; name, systems)
end
