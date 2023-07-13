"""
    Fixed(;name, s_0=0.0)

Flange fixed in housing at a given position.

# Parameters:

  - `s_0`: [m] Fixed offset position of housing

# Connectors:

  - `flange: 1-dim. translational flange`
"""
@mtkmodel Fixed begin#(; name, s_0 = 0.0)
    @parameters begin
        s_0
    end
    @components begin
        flange = Flange(; s = s_0)
    end
    @equations begin
        flange.s ~ s_0
    end
end

"""
    Mass(; name, m, s = 0.0, v = 0.0)

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
@mtkmodel Mass begin
    @parameters begin
        m
    end
    @variables begin
        s(t) = 0.0
        v(t) = 0.0
        f(t) = 0.0
    end
    @components begin
        flange = Flange(; s = s)
    end
    @equations begin
        flange.s ~ s
        flange.f ~ f
        D(s) ~ v
        D(v) ~ f / m
    end
end

const REL = Val(:relative)
@component function Spring(::Val{:relative}; name, k, v_a_0 = 0.0, v_b_0 = 0.0,
    delta_s_0 = 0,
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

@component function Spring(::Val{:absolute}; name, k, s_a_0 = 0, s_b_0 = 0, l = 0)
    pars = @parameters begin
        k = k
        s_a_0 = s_a_0
        s_b_0 = s_b_0
        l = l
    end
    vars = @variables begin
        # delta_s(t) = s1 - s2
        f(t) = k * (s_a_0 - s_b_0 - l)
    end

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
    Damper(; name, d, v1 =0.0, v2 = 0.0, flange_a.s = 0, flange_b.s = 0)

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
@mtkmodel Damper begin#(; name, d, v_a_0 = 0.0, v_b_0 = 0.0, s_a_0 = 0, s_b_0 = 0)
    @parameters begin
        d
    end
    @variables begin
        v1(t) = 0.0
        v2(t) = 0.0
        f(t) = +(v1 - v2) * d
    end

    @components begin
        flange_a = Flange(; s = 0.0, f = (v1 - v2) * d)
        flange_b = Flange(; s = 0.0, f = -(v1 - v2) * d)
    end

    @equations begin
        D(flange_a.s) ~ v1
        D(flange_b.s) ~ v2
        f ~ (v1 - v2) * d
        flange_a.f ~ +f
        flange_b.f ~ -f
    end
end
