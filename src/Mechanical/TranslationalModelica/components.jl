"""
    Fixed(;name, s0=0.0)

Flange fixed in housing at a given position.

# Parameters:

  - `s0`: [m] Fixed offset position of housing

# Connectors:

  - `flange: 1-dim. translational flange`
"""
function Fixed(; name, s0 = 0.0)
    pars = @parameters s0 = s0
    vars = []

    @named flange = Flange()

    eqs = [flange.s ~ s0]

    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [flange.s => s0]),
                   flange)
end

"""
    Mass(; name, m, s0 = 0.0, v0 = 0.0)

Sliding mass with inertia

# Parameters:

  - `m`: [kg] Mass of sliding mass
  - `s0`: [m] Initial value of absolute position of sliding mass
  - `v0`: [m/s] Initial value of absolute linear velocity of sliding mass

# States:

  - `s`: [m] Absolute position of sliding mass
  - `v`: [m/s] Absolute linear velocity of sliding mass (= D(s))

# Connectors:

  - `flange: 1-dim. translational flange of mass`
"""
function Mass(m; name, s0 = 0.0, v0 = 0.0)
    @named pr = PartialRigid(; L=0, s0)
    @unpack flange_a, flange_b, s = pr
    @parameters m = m [description = "Mass of sliding mass [kg]"]
    @variables v(t)=v0 [description = "Absolute linear velocity of sliding mass [m/s]"]
    @variables a(t)=0 [description = "Absolute linear acceleration of sliding mass [m/s^2]"]
    eqs = [
        v ~ D(s)
        a ~ D(v)
        m*a ~ flange_a.f + flange_b.f
    ]
    return extend(ODESystem(eqs, t; name), pr)
end


"""
    Spring(c; name, s_rel0=0)

Linear 1D translational spring

# Parameters:

  - `c`: [N/m] Spring constant
  - `s_rel0`: Unstretched spring length

# Connectors:

  - `flange_a: 1-dim. translational flange on one side of spring`
  - `flange_b: 1-dim. translational flange on opposite side of spring` #default function
"""
function Spring(c; name, s_rel0 = 0)
    @named pc = PartialCompliant()
    @unpack flange_a, flange_b, s_rel, f = pc
    @parameters c = c [description = "Spring constant [N/m]"]
    @parameters s_rel0 = s_rel0 [description = "Unstretched spring length [m]"]

    eqs = [f ~ c*(s_rel - s_rel0)]
    return extend(ODESystem(eqs, t; name), pc)
end

"""
    Damper(d; name)

Linear 1D translational damper

# Parameters:

  - `d`: [N.s/m] Damping constant

# Connectors:

  - `flange_a: 1-dim. translational flange on one side of damper`
  - `flange_b: 1-dim. translational flange on opposite side of damper`
"""
function Damper(d; name)
  @named pc = PartialCompliantWithRelativeStates()
  @unpack flange_a, flange_b, v_rel, f = pc
  @parameters d = d [description = "Damping constant [Ns/m]"]
  @variables lossPower(t) = 0 [description = "Power dissipated by the damper [W]"]
  eqs = [f ~ d*v_rel; lossPower ~ f*v_rel]
  return extend(ODESystem(eqs, t; name), pc)
end
