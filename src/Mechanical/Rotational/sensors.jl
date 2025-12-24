"""
    AngleSensor(; name)

Ideal sensor to measure the absolute flange angle

# Connectors:

  - `flange`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `phi`: [RealOutput](@ref) Absolute angle of flange
"""
@component function AngleSensor(; name)
    pars = @parameters begin
    end

    systems = @named begin
        flange = Flange()
        phi = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        phi.u ~ flange.phi,
        flange.tau ~ 0
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    SpeedSensor(; name)

Ideal sensor to measure the absolute flange angular velocity

# Connectors:

  - `flange`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `w`: [RealOutput](@ref) Absolute angular velocity of flange
"""
@component function SpeedSensor(; name)
    pars = @parameters begin
    end

    systems = @named begin
        flange = Flange()
        w = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        D(flange.phi) ~ w.u,
        flange.tau ~ 0
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    TorqueSensor(;name)

Ideal sensor to measure the torque between two flanges (`= flange_a.tau`)

# Connectors:

  - `flange_a`: [Flange](@ref) Left flange of shaft
  - `flange_b`: [Flange](@ref) Left flange of shaft
  - `tau`: [RealOutput](@ref) Torque in flange flange_a and flange_b (`tau = flange_a.tau = -flange_b.tau`)
"""
@component function TorqueSensor(; name)
    pars = @parameters begin
    end

    systems = @named begin
        flange_a = Flange()
        flange_b = Flange()
        tau = RealOutput()
    end

    vars = @variables begin
    end

    equations = Equation[
        flange_a.phi ~ flange_b.phi,
        tau.u ~ flange_a.tau
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    RelSpeedSensor(; name)

Ideal sensor to measure the relative angular velocity

# Connectors:

  - `flange_a`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `flange_b`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `w`: [RealOutput](@ref) Absolute angular velocity of flange
"""
@component function RelSpeedSensor(; name, phi_rel = nothing)
    pars = @parameters begin
    end

    systems = @named begin
        flange_a = Flange()
        flange_b = Flange()
        w_rel = RealOutput()
    end

    vars = @variables begin
        phi_rel(t) = phi_rel, [guess = 0.0]
    end

    equations = Equation[
        0 ~ flange_a.tau + flange_b.tau,
        phi_rel ~ flange_b.phi - flange_a.phi,
        D(phi_rel) ~ w_rel.u,
        0 ~ flange_a.tau
    ]

    return System(equations, t, vars, pars; name, systems)
end
