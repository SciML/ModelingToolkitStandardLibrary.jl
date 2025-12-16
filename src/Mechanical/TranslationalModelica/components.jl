"""
    Fixed(; name, s0 = 0.0)

Flange fixed in housing at a given position.

# Parameters:

  - `s0`: [m] Fixed offset position of housing

# Connectors:

  - `flange: 1-dim. translational flange`
"""
@component function Fixed(; name, s0 = 0)
    pars = @parameters begin
        s0 = s0
    end

    systems = @named begin
        flange = Flange()
    end

    vars = @variables begin
    end

    equations = Equation[
        flange.s ~ s0
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Mass(; name, m, s, v = 0.0)

Sliding mass with inertia

# Parameters:

  - `m`: [kg] Mass of sliding mass

# States:

  - `s`: [m] Absolute position of sliding mass. It accepts an initial value, which defaults to 0.0.
  - `v`: [m/s] Absolute linear velocity of sliding mass (= D(s)). It accepts an initial value, which defaults to 0.0.

# Connectors:

  - `flange: 1-dim. translational flange of mass`
"""
@component function Mass(; name, m = 0.0, s = nothing, v = nothing, a = nothing)
    @named pr = PartialRigid(; L = 0.0, s)
    @unpack flange_a, flange_b, s = pr

    pars = @parameters begin
        m = m, [description = "Mass of sliding mass [kg]"]
    end

    systems = @named begin
    end

    vars = @variables begin
        v(t) = v, [description = "Absolute linear velocity of sliding mass [m/s]"]
        a(t) = a, [description = "Absolute linear acceleration of sliding mass [m/s^2]"]
    end

    equations = Equation[
        v ~ D(s),
        a ~ D(v),
        m * a ~ flange_a.f + flange_b.f
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, pr)
end

"""
    Spring(; c= 0.0, name, s_rel0 = 0)

Linear 1D translational spring

# Parameters:

  - `c`: [N/m] Spring constant
  - `s_rel0`: Unstretched spring length

# Connectors:

  - `flange_a: 1-dim. translational flange on one side of spring`
  - `flange_b: 1-dim. translational flange on opposite side of spring` #default function
"""
@component function Spring(; name, c = 0.0, s_rel0 = 0.0)
    @named pc = PartialCompliant()
    @unpack flange_a, flange_b, s_rel, f = pc

    pars = @parameters begin
        c = c, [description = "Spring constant [N/m]"]
        s_rel0 = s_rel0, [description = "Unstretched spring length [m]"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        f ~ c * (s_rel - s_rel0)
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, pc)
end

"""
    Damper(; name, d = 0.0)

Linear 1D translational damper

# Parameters:

  - `d`: [N.s/m] Damping constant

# Connectors:

  - `flange_a: 1-dim. translational flange on one side of damper`
  - `flange_b: 1-dim. translational flange on opposite side of damper`
"""
@component function Damper(; name, d = 0.0, lossPower = nothing)
    @named pc = PartialCompliantWithRelativeStates()
    @unpack flange_a, flange_b, v_rel, f = pc

    pars = @parameters begin
        d = d, [description = "Damping constant [Ns/m]"]
    end

    systems = @named begin
    end

    vars = @variables begin
        lossPower(t) = lossPower, [description = "Power dissipated by the damper [W]"]
    end

    equations = Equation[
        f ~ d * v_rel,
        lossPower ~ f * v_rel
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, pc)
end

"""
    SpringDamper(; name, c = 0.0, d = 0.0, s_rel0 = 0.0)

Linear 1D translational spring and damper in parallel

# Parameters:
- `c`: [N/m] Spring constant
- `d`: [N.s/m] Damping constant
- `s_rel0`: Unstretched spring length

# Connectors:
- `flange_a: 1-dim. translational flange on one side of spring`
- `flange_b: 1-dim. translational flange on opposite side of spring`

# Variables:
- `lossPower`: [W] Power dissipated by the damper
- `f`: [N] Total force
"""
@component function SpringDamper(; name, c = 0.0, d = 0.0, s_rel0 = 0.0, lossPower = nothing)
    @named pc = PartialCompliantWithRelativeStates()
    @unpack flange_a, flange_b, s_rel, v_rel, f = pc

    pars = @parameters begin
        d = d, [description = "Damping constant [Ns/m]"]
        c = c, [description = "Spring constant [N/m]"]
        s_rel0 = s_rel0, [description = "Unstretched spring length [m]"]
    end

    systems = @named begin
    end

    vars = @variables begin
        lossPower(t) = lossPower, [description = "Power dissipated by the damper [W]"]
    end

    equations = Equation[
        f ~ c * (s_rel - s_rel0) + d * v_rel,
        lossPower ~ d * v_rel^2
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, pc)
end
