"""
    Free(; name)

Use to close a system that has un-connected `MechanicalPort`'s where the force should not be zero (i.e. you want to solve for the force to produce the given movement of the port)

# Connectors:

  - `flange`: 1-dim. translational flange
"""
@mtkmodel Free begin
    @components begin
        flange = MechanicalPort()
    end
    @variables begin
        f(t) = 0.0
    end
    @equations begin
        flange.f ~ -f
    end
end

"""
    Fixed(; name)

Fixes a flange position (velocity = 0)

# Connectors:

  - `flange`: 1-dim. translational flange
"""
@mtkmodel Fixed begin
    @components begin
        flange = MechanicalPort()
    end
    @equations begin
        flange.v ~ 0
    end
end

"""
    Mass(; name, v_0 = 0.0, m, s = nothing, g = nothing)

Sliding mass with inertia

# Parameters:

  - `m`: [kg] mass of sliding body
  - `v_0`: [m/s] Initial value of absolute linear velocity of sliding mass (default 0 m/s)
  - `s`: [m] (optional) initial value of absolute position of sliding mass
  - `g`: [m/s²] (optional) gravity field acting on the mass, positive value acts in the positive direction

# States:

  - `v`: [m/s] absolute linear velocity of sliding mass
  - `s`: [m] (optional with parameter s) absolute position of sliding mass

# Connectors:

  - `flange`: 1-dim. translational flange
"""
@component function Mass(; name, v = 0.0, m, s = nothing, g = nothing)
    pars = @parameters begin
        m = m
    end
    @named flange = MechanicalPort(; v = v)

    vars = @variables begin
        v(t) = v
        f(t) = 0
    end

    eqs = [flange.v ~ v
        flange.f ~ -f]

    # gravity option
    if g !== nothing
        @parameters g = g
        push!(pars, g)
        push!(eqs, D(v) ~ f / m + g)
    else
        push!(eqs, D(v) ~ f / m)
    end

    # position option
    if s !== nothing
        @variables s(t) = s
        push!(vars, s)

        push!(eqs, D(s) ~ v)
    end

    return compose(ODESystem(eqs, t, vars, pars; name = name),
        flange)
end

const REL = Val(:relative)

"""
    Spring(; name, k, delta_s = 0.0,  va=0.0, v_b_0=0.0)

Linear 1D translational spring

# Parameters:

  - `k`: [N/m] Spring constant
  - `delta_s`: initial spring stretch
  - `va`: [m/s] Initial value of absolute linear velocity at flange_a (default 0 m/s)
  - `v_b_0`: [m/s] Initial value of absolute linear velocity at flange_b (default 0 m/s)

# Connectors:

  - `flange_a`: 1-dim. translational flange on one side of spring
  - `flange_b`: 1-dim. translational flange on opposite side of spring
"""
@component function Spring(; name, k, delta_s = 0.0, flange_a__v = 0.0, flange_b__v = 0.0)
    Spring(REL; name, k, delta_s, flange_a__v, flange_b__v)
end # default

@component function Spring(::Val{:relative}; name, k, delta_s = 0.0, flange_a__v = 0.0,
        flange_b__v = 0.0)
    pars = @parameters begin
        k = k
    end
    vars = @variables begin
        delta_s(t) = delta_s
        f(t) = 0
    end

    @named flange_a = MechanicalPort(; v = flange_a__v)
    @named flange_b = MechanicalPort(; v = flange_b__v)

    eqs = [D(delta_s) ~ flange_a.v - flange_b.v
        f ~ k * delta_s
        flange_a.f ~ -f
        flange_b.f ~ +f]
    return compose(ODESystem(eqs, t, vars, pars; name = name),
        flange_a,
        flange_b) #flange_a.f => +k*delta_s, flange_b.f => -k*delta_s
end

const ABS = Val(:absolute)
@component function Spring(::Val{:absolute}; name, k, sa = 0, sb = 0, flange_a__v = 0.0,
        flange_b__v = 0.0, l = 0)
    pars = @parameters begin
        k = k
        l = l
    end
    vars = @variables begin
        sa(t) = sa
        sb(t) = sb
        f(t) = 0
    end

    @named flange_a = MechanicalPort(; v = flange_a__v)
    @named flange_b = MechanicalPort(; v = flange_b__v)

    eqs = [D(sa) ~ flange_a.v
        D(sb) ~ flange_b.v
        f ~ k * (sa - sb - l) #delta_s
        flange_a.f ~ -f
        flange_b.f ~ +f]
    return compose(ODESystem(eqs, t, vars, pars; name = name),
        flange_a,
        flange_b) #, flange_a.f => k * (flange_a__s - flange_b__s - l)
end

"""
    Damper(; name, d, flange_a.v = 0.0, flange_b.v = 0.0)

Linear 1D translational damper

# Parameters:

  - `d`: [N.s/m] Damping constant

# Connectors:

  - `flange_a`: 1-dim. translational flange on one side of damper. Initial value of state `v` is set to 0.0 m/s.
  - `flange_b`: 1-dim. translational flange on opposite side of damper. Initial value of state `v` is set to 0.0 m/s.
"""
@mtkmodel Damper begin
    @parameters begin
        d
    end
    @variables begin
        v(t)
        f(t) = 0.0
    end

    @components begin
        flange_a = MechanicalPort(; v = 0.0)
        flange_b = MechanicalPort(; v = 0.0)
    end

    @equations begin
        v ~ flange_a.v - flange_b.v
        f ~ v * d
        flange_a.f ~ -f
        flange_b.f ~ +f
    end
end
