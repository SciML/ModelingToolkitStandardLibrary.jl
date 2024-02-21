"""
    Force(; name)

Linear 1D force input source

# Connectors:

  - `flange`: 1-dim. translational flange
  - `f`: real input
"""
@mtkmodel Force begin
    @components begin
        flange = MechanicalPort(; v = 0.0)
        f = RealInput(unit = u"N")
    end

    @equations begin
        flange.f ~ -f.u
    end
end

"""
    Position(solves_force = true; name)

Linear 1D position input source.  Set `solves_force=false` to force input force to 0 (i.e. only the position is given, the respective force needed is already provided elsewhere in the model).

# Connectors:

  - `flange`: 1-dim. translational flange
  - `s`: real input
"""
@mtkmodel Position begin
    @structural_parameters begin
        solves_force = true
    end

    @components begin
        flange = MechanicalPort(; v = 0)
        s = RealInput(unit = u"m")
    end

    @equations begin
        D(s.u) ~ flange.v
        0 ~ flange.f
    end
end

"""
    Velocity(solves_force = true; name)

Linear 1D position input source.  Set `solves_force=false` to force input force to 0 (i.e. only the velocity is given, the respective force needed is already provided elsewhere in the model).

# Connectors:

  - `flange`: 1-dim. translational flange
  - `v`: real input
"""
@component function Velocity(solves_force = true; name)
    systems = @named begin
        flange = MechanicalPort(; v = 0)
        v = RealInput(unit = u"m/s")
    end

    eqs = [
        v.u ~ flange.v
    ]

    !solves_force && push!(eqs, 0 ~ flange.f)

    ODESystem(eqs, t, [], []; name, systems)
end

"""
Acceleration(solves_force = true; name)

Linear 1D position input source.  Set `solves_force=false` to force input force to 0 (i.e. only the acceleration is given, the respective force needed is already provided elsewhere in the model).

# Connectors:

  - `flange`: 1-dim. translational flange
  - `a`: real input
"""

@mtkmodel Acceleration begin#(solves_force = true; name)
    @structural_parameters begin
        solves_force = true
    end
    @variables begin
        v(t) = 0, [unit = u"m/s"]
    end
    @components begin
        flange = MechanicalPort(; v = 0)
        a = RealInput(unit = u"m/s^2")
    end
    @equations begin
        v ~ flange.v
        D(v) ~ a.u
        if solves_force
            0 ~ flange.f
        end
    end
end
