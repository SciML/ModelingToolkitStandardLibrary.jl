"""
    Fixed(;name, s_0=0.0)

Flange fixed in housing at a given position.

# Parameters:

  - `s_0`: [m] Fixed offset position of housing

# Connectors:

  - `flange: 1-dim. translational flange`
"""
@mtkmodel Fixed begin
    @parameters begin
        s_0, [description = "Fixed offset position of housing", unit = u"m"]
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

# States:

  - `s`: [m] Absolute position of sliding mass. It accepts an initial value, which defaults to 0.0.
  - `v`: [m/s] Absolute linear velocity of sliding mass (= der(s)). It accepts an initial value, which defaults to 0.0.
  - `f`: [N] Force. It accepts an initial value, which defaults to 0.0.

# Connectors:

  - `flange: 1-dim. translational flange of mass`
"""
@mtkmodel Mass begin
    @parameters begin
        m, [description = "Mass of sliding mass", unit = u"kg"]
    end
    @variables begin
        s(t) = 0.0, [description = "Absolute position of sliding mass", unit = u"m"]
        function v(t)
            0.0,
            [description = "Absolute linear velocity of sliding mass", unit = u"m*s^-1"]
        end
        f(t) = 0.0, [description = "Force", unit = u"N"]
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
@component function Spring(::Val{:relative}; name, k, va = 0.0, vb = 0.0,
    delta_s = 0, flange_a__s = 0, flange_b__s = 0)
    pars = @parameters begin
        k = k
    end
    vars = @variables begin
        va(t) = va
        vb(t) = vb
        delta_s(t) = delta_s
        f(t) = delta_s * k
    end

    @named flange_a = Flange()
    @named flange_b = Flange()

    eqs = [D(flange_a.s) ~ va
        D(flange_b.s) ~ vb
        D(delta_s) ~ va - vb
        f ~ k * delta_s
        flange_a.f ~ +f
        flange_b.f ~ -f]

    return compose(ODESystem(eqs, t, vars, pars; name = name,
            defaults = [
                flange_a.s => flange_a__s,
                flange_b.s => flange_b__s,
                flange_a.f => +delta_s * k,
                flange_b.f => -delta_s * k,
            ]), flange_a, flange_b)
end

const ABS = Val(:absolute)

"""
    Spring(; name, k, flange_a__s = 0, flange_b__s = 0, l=0)

Linear 1D translational spring

# Parameters:

  - `k`: [N/m] Spring constant
  - `l`: Unstretched spring length
  - `flange_a__s`: [m] Initial value of absolute position of flange_a
  - `flange_b__s`: [m] Initial value of absolute position of flange_b

# Connectors:

  - `flange_a: 1-dim. translational flange on one side of spring`
  - `flange_b: 1-dim. translational flange on opposite side of spring` #default function
"""
function Spring(; name, k, flange_a__s = 0, flange_b__s = 0, l = 0)
    Spring(ABS; name, k, flange_a__s, flange_b__s, l)
end #default function

@component function Spring(::Val{:absolute};
    name, k, flange_a__s = 0,
    flange_b__s = 0, l = 0)
    pars = @parameters begin
        k = k
        l = l
    end
    vars = @variables begin
        f(t) = k * (flange_a__s - flange_b__s - l)
    end

    @named flange_a = Flange(; s = flange_a__s, f = k * (flange_a__s - flange_b__s - l))
    @named flange_b = Flange(; s = flange_a__s, f = -k * (flange_a__s - flange_b__s - l))

    eqs = [
    #    delta_s ~ flange_a.s - flange_b.s
        f ~ k * (flange_a.s - flange_b.s - l) #delta_s
        flange_a.f ~ +f
        flange_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name), flange_a, flange_b)
end

"""
    Damper(; name, d, va =0.0, vb = 0.0, flange_a.s = 0, flange_b.s = 0)

Linear 1D translational damper

# Parameters:

  - `d`: [N.s/m] Damping constant
  - `flange_a__s`: [m] Initial value of absolute position of flange_a
  - `flange_b__s`: [m] Initial value of absolute position of flange_b

# Connectors:

  - `flange_a: 1-dim. translational flange on one side of damper`
  - `flange_b: 1-dim. translational flange on opposite side of damper`
"""
@mtkmodel Damper begin
    @parameters begin
        d, [description = "Damping constant", unit = u"N.s/m"]
    end
    @variables begin
        va(t) = 0.0
        vb(t) = 0.0
        f(t) = +(va - vb) * d
    end

    @components begin
        flange_a = Flange(; s = 0.0, f = (va - vb) * d)
        flange_b = Flange(; s = 0.0, f = -(va - vb) * d)
    end

    @equations begin
        D(flange_a.s) ~ va
        D(flange_b.s) ~ vb
        f ~ (va - vb) * d
        flange_a.f ~ +f
        flange_b.f ~ -f
    end
end
