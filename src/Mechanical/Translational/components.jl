"""
    Fixed(;name, s0=0.0)

Flange fixed in housing at a given position.

# Parameters:
- `s0`: [m] Fixed offset position of housing

# Connectors:
- `flange: 1-dim. translational flange`
"""
function Fixed(; name, s0 = 0.0)
    @named flange = Flange()
    @parameters s0 = s0
    eqs = [flange.s ~ s0]
    return compose(ODESystem(eqs, t, [], [s0]; name = name), flange)
end

"""
    Mass(;name, m, s_start=0.0, v_start=0.0, a_start=0.0)

Sliding mass with inertia

# Parameters:
- `m`: [kg] Mass of sliding mass
- `s_start`: [m] Initial value of absolute position of sliding mass
- `v_start`: [m/s] Initial value of absolute linear velocity of sliding mass
- `a_start`: [m/s²] Initial value of absolute linear acceleration of sliding mass

# States: 
- `s`: [m] Absolute position of sliding mass
- `v`: [m/s] Absolute linear velocity of sliding mass (= der(s)) 
- `a`: [m/s²] Absolute linear acceleration of sliding mass (= der(v))

# Connectors:
- `flange_a: 1-dim. translational flange on one side of mass`
- `flange_b: 1-dim. translational flange on opposite side of mass`
"""
function Mass(; name, m, s_start = 0.0, v_start = 0.0, a_start = 0.0)
    @named flange_a = Flange()
    @named flange_b = Flange()
    @parameters m = m
    sts = @variables begin
        s(t) = s_start
        v(t) = v_start
        a(t) = a_start
    end
    eqs = [s ~ flange_a.s
           s ~ flange_b.s
           D(s) ~ v
           D(v) ~ a
           m * a ~ flange_a.f + flange_b.f]
    return compose(ODESystem(eqs, t, sts, [m]; name = name), flange_a, flange_b)
end

"""
    Spring(;name, c, s_rel0=0.0)

Linear 1D translational spring

# Parameters:
- `c`: [N/m] Spring constant
- `s_rel0`: Unstretched spring length

# Connectors:
- `flange_a: 1-dim. translational flange on one side of spring`
- `flange_b: 1-dim. translational flange on opposite side of spring`
"""
function Spring(; name, c, s_rel0 = 0.0)
    @named partial_comp = PartialCompliant()
    @unpack s_rel, f = partial_comp
    pars = @parameters begin
        c = c
        s_rel0 = s_rel0
    end
    eqs = [f ~ c * (s_rel - s_rel0)]
    extend(ODESystem(eqs, t, [], pars; name = name), partial_comp)
end

"""
    Damper(;name, d) 

Linear 1D translational damper

# Parameters:
- `d`: [N.s/m] Damping constant

# Connectors:
- `flange_a: 1-dim. translational flange on one side of damper`
- `flange_b: 1-dim. translational flange on opposite side of damper`
"""
function Damper(; name, d)
    @named partial_comp = PartialCompliantWithRelativeStates()
    @unpack v_rel, f = partial_comp
    pars = @parameters d = d
    eqs = [f ~ d * v_rel]
    extend(ODESystem(eqs, t, [], pars; name = name), partial_comp)
end
