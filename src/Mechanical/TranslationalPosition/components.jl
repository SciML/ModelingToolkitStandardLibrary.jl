"""
    Fixed(;name, s_0=0.0)

Flange fixed in housing at a given position.

# Parameters:

  - `s_0`: [m] Fixed offset position of housing

# Connectors:

  - `flange: 1-dim. translational flange`
"""
function Fixed(; name, s_0 = 0.0)
    pars = @parameters s_0 = s_0
    vars = []

    @named flange = Flange()

    eqs = [flange.s ~ s_0]

    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [flange.s => s_0]),
                   flange)
end

"""
    Mass(; name, m, s_0 = 0.0, v_0 = 0.0)

Sliding mass with inertia

# Parameters:

  - `m`: [kg] Mass of sliding mass
  - `s_0`: [m] Initial value of absolute position of sliding mass
  - `v_0`: [m/s] Initial value of absolute linear velocity of sliding mass

# States:

  - `s`: [m] Absolute position of sliding mass
  - `v`: [m/s] Absolute linear velocity of sliding mass (= der(s))

# Connectors:

  - `flange: 1-dim. translational flange of mass`
"""
function Mass(; name, m, s_0 = 0.0, v_0 = 0.0)
    @named flange = Flange()
    pars = @parameters begin
        m = m
        s_0 = s_0
        v_0 = v_0
    end
    vars = @variables begin
        s(t) = s_0
        v(t) = v_0
        f(t) = 0
    end
    eqs = [flange.s ~ s
           flange.f ~ f
           D(s) ~ v
           D(v) ~ f / m]
    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [flange.s => s_0]),
                   flange)
end

const REL = Val(:relative)
function Spring(::Val{:relative}; name, k, v_a_0 = 0.0, v_b_0 = 0.0, delta_s_0 = 0,
                s_a_0 = 0,
                s_b_0 = 0)
    pars = @parameters begin
        k = k
        v_a_0 = v_a_0
        v_b_0 = v_b_0
        s_a_0 = s_a_0
        s_b_0 = s_b_0
        delta_s_0 = delta_s_0
    end
    vars = @variables begin
        v1(t) = v_a_0
        v2(t) = v_b_0
        delta_s(t) = delta_s_0
        f(t) = delta_s_0 * k
    end

    @named flange_a = Flange()
    @named flange_b = Flange()

    eqs = [D(flange_a.s) ~ v1
           D(flange_b.s) ~ v2
           D(delta_s) ~ v1 - v2
           f ~ k * delta_s
           flange_a.f ~ +f
           flange_b.f ~ -f]

    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [
                                 flange_a.s => s_a_0,
                                 flange_b.s => s_b_0,
                                 flange_a.f => +delta_s_0 * k,
                                 flange_b.f => -delta_s_0 * k,
                             ]), flange_a, flange_b)
end

const ABS = Val(:absolute)

"""
    Spring(; name, k, s_a_0 = 0, s_b_0 = 0, l=0)

Linear 1D translational spring

# Parameters:

  - `k`: [N/m] Spring constant
  - `l`: Unstretched spring length
  - `s_a_0`: [m] Initial value of absolute position of flange_a
  - `s_b_0`: [m] Initial value of absolute position of flange_b

# Connectors:

  - `flange_a: 1-dim. translational flange on one side of spring`
  - `flange_b: 1-dim. translational flange on opposite side of spring` #default function
"""
Spring(; name, k, s_a_0 = 0, s_b_0 = 0, l = 0) = Spring(ABS; name, k, s_a_0, s_b_0, l) #default function

function Spring(::Val{:absolute}; name, k, s_a_0 = 0, s_b_0 = 0, l = 0)
    pars = @parameters begin
        k = k
        s_a_0 = s_a_0
        s_b_0 = s_b_0
        l = l
    end
    vars = @variables begin
    # delta_s(t) = s1 - s2
    f(t) = k * (s_a_0 - s_b_0 - l) end

    @named flange_a = Flange()
    @named flange_b = Flange()

    eqs = [
           #    delta_s ~ flange_a.s - flange_b.s
           f ~ k * (flange_a.s - flange_b.s - l) #delta_s
           flange_a.f ~ +f
           flange_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [
                                 flange_a.s => s_a_0,
                                 flange_b.s => s_b_0,
                                 flange_a.f => k * (s_a_0 - s_b_0 - l),
                                 flange_b.f => -k * (s_a_0 - s_b_0 - l),
                             ]), flange_a, flange_b)
end

"""
    Damper(; name, d, v_a_0=0.0, v_b_0=0.0, s_a_0 = 0, s_b_0 = 0)

Linear 1D translational damper

# Parameters:

  - `d`: [N.s/m] Damping constant
  - `s_a_0`: [m] Initial value of absolute position of flange_a
  - `s_b_0`: [m] Initial value of absolute position of flange_b
  - `v_a_0`: [m/s] Initial value of absolute linear velocity of flange_a
  - `v_b_0`: [m/s] Initial value of absolute linear velocity of flange_a

# Connectors:

  - `flange_a: 1-dim. translational flange on one side of damper`
  - `flange_b: 1-dim. translational flange on opposite side of damper`
"""
function Damper(; name, d, v_a_0 = 0.0, v_b_0 = 0.0, s_a_0 = 0, s_b_0 = 0)
    pars = @parameters begin
        d = d
        s_a_0 = s_a_0
        s_b_0 = s_b_0
        v_a_0 = v_a_0
        v_b_0 = v_b_0
    end
    vars = @variables begin
        v1(t) = v_a_0
        v2(t) = v_b_0
        f(t) = +(v_a_0 - v_b_0) * d
    end

    @named flange_a = Flange()
    @named flange_b = Flange()

    eqs = [D(flange_a.s) ~ v1
           D(flange_b.s) ~ v2
           f ~ (v1 - v2) * d
           flange_a.f ~ +f
           flange_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [
                                 flange_a.s => s_a_0,
                                 flange_b.s => s_b_0,
                                 flange_a.f => +(v_a_0 - v_b_0) * d,
                                 flange_b.f => -(v_a_0 - v_b_0) * d,
                             ]), flange_a, flange_b)
end
