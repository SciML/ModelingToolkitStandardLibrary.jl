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
        f = RealInput()
    end

    @equations begin
        flange.f ~ -f.u
    end
end

"""
    Position(; s.u_start = 0.0, name)

Linear 1D position input source

# Connectors:

  - `flange`: 1-dim. translational flange
  - `s`: real input. `s.u_start` accepts initial value and defaults to 0.0.
"""
@mtkmodel Position begin
    @parameters begin
        solves_force = false
    end
    @components begin
        flange = MechanicalPort(; v = 0.0)
        s = RealInput(; u_start = 0.0)
    end
    @variables begin
        x(t)
    end
    @equations begin
        D(x) ~ flange.v
        s.u ~ x
        !getdefault(solves_force) && 0 ~ flange.f
    end
end

@mtkmodel Velocity begin
    @parameters begin
        solves_force = false
    end
    @components begin
        flange = MechanicalPort()
        v = RealInput()
    end
    @equations begin
        v.u ~ flange.v
        getdefault(solves_force) && 0 ~ flange.f
    end
end

@mtkmodel Acceleration begin
    @parameters begin
        solves_force = false
    end
    @components begin
        flange = MechanicalPort(; v = 0.0)
        a = RealInput()
    end
    @variables begin
        v(t) = 0
    end

    @equations begin
        v ~ flange.v
        D(v) ~ a.u
        !getdefault(solves_force) && 0 ~ flange.f
    end
end
