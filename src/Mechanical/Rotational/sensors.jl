"""
    AngleSensor(;name)

Ideal sensor to measure the absolute flange angle

# Connectors:

  - `flange`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `phi`: [RealOutput](@ref) Absolute angle of flange
"""
@component function AngleSensor(; name)
    @named flange = Flange()
    @named phi = RealOutput()
    eqs = [phi.u ~ flange.phi
           flange.tau ~ 0]
    return ODESystem(eqs, t, [], []; name = name, systems = [flange, phi])
end

"""
    SpeedSensor(;name)

Ideal sensor to measure the absolute flange angular velocity

# Connectors:

  - `flange`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `w`: [RealOutput](@ref) Absolute angular velocity of flange
"""
@component function SpeedSensor(; name)
    @named flange = Flange()
    @named w = RealOutput()
    eqs = [D(flange.phi) ~ w.u
           flange.tau ~ 0]
    return ODESystem(eqs, t, [], []; name = name, systems = [flange, w])
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
    @named flange_a = Flange()
    @named flange_b = Flange()
    @named tau = RealOutput()
    eqs = [flange_a.phi ~ flange_b.phi
           tau.u ~ flange_a.tau]
    return ODESystem(eqs, t, [], []; name = name, systems = [flange_a, flange_b, tau])
end

"""
    RelSpeedSensor(;name)

Ideal sensor to measure the relative angular velocity

# Connectors:

  - `flange_a`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `flange_b`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `w`: [RealOutput](@ref) Absolute angular velocity of flange
"""
@component function RelSpeedSensor(; name)
    @named flange_a = Flange()
    @named flange_b = Flange()
    @named w_rel = RealOutput()
    @variables phi_rel(t) = 0.0
    eqs = [0 ~ flange_a.tau + flange_b.tau
           phi_rel ~ flange_b.phi - flange_a.phi
           D(phi_rel) ~ w_rel.u
           0 ~ flange_a.tau]
    return ODESystem(eqs, t, [phi_rel], []; name = name,
                     systems = [flange_a, flange_b, w_rel])
end
