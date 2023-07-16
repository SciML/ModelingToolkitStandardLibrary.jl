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
  - `s`: real input. `s.u_start` accepts an initial value, which It accepts an initial value, which defaults to 0.0.
"""
@component function Position(; solves_force = true, s__u_start = 0, name)
    systems = @named begin
        flange = MechanicalPort()
        s = RealInput()
    end

    vars = @variables x(t)

    eqs = [D(x) ~ flange.v
        s.u ~ x]

    !solves_force && push!(eqs, 0 ~ flange.f)

    ODESystem(eqs, t, vars, [];
        name, systems, defaults = [flange.v => 0, s.u => s__u_start])
end

@component function Velocity(; solves_force = true, name)
    systems = @named begin
        flange = MechanicalPort()
        v = RealInput()
    end

    eqs = [
        v.u ~ flange.v,
    ]

    !solves_force && push!(eqs, 0 ~ flange.f)

    ODESystem(eqs, t, [], []; name, systems, defaults = [flange.v => 0])
end

@component function Acceleration(solves_force = true; s__u_start = 0, name)
    systems = @named begin
        flange = MechanicalPort()
        a = RealInput(; u_start = s__u_start)
    end

    vars = @variables v(t) = 0

    eqs = [v ~ flange.v
        D(v) ~ a.u]

    !solves_force && push!(eqs, 0 ~ flange.f)

    ODESystem(eqs, t, vars, []; name, systems, defaults = [flange.v => 0])
end
