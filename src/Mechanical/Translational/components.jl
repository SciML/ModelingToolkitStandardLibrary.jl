"""
    Fixed(;name)

Fixes a flange position (velocity = 0)

# Connectors:
- `flange: 1-dim. translational flange`
"""
function Fixed(; name)
    @named flange = MechanicalPort()
    eqs = [flange.v ~ 0]
    return compose(ODESystem(eqs, t, [], []; name = name, defaults = [flange.v => 0]), flange)
end

"""
    Mass(; name, v_0 = 0.0, m, s_0 = nothing, g = nothing)

Sliding mass with inertia

# Parameters:
- `m`: [kg] mass of sliding body
- `v_0`: [m/s] Initial value of absolute linear velocity of sliding mass (default 0 m/s)
- `s_0`: [m] (optional) initial value of absolute position of sliding mass
- `g`: [m/s²] (optional) gravity field acting on the mass, positive value acts in the positive direction


# States:
- `v`: [m/s] absolute linear velocity of sliding mass
- `s`: [m] (optional with parameter s_0) absolute position of sliding mass


# Connectors:
- `flange: 1-dim. translational flange`
"""
function Mass(; name, v_0 = 0.0, m, s_0 = nothing, g = nothing)
    pars = @parameters begin
        m = m
        v_0 = v_0
    end
    vars = @variables begin
        v(t) = v_0
        f(t) = 0
    end

    @named flange = MechanicalPort()

    eqs = [flange.v ~ v
           flange.f ~ f]

    # gravity option
    if !isnothing(g)
        @parameters g = g
        push!(pars, g)
        push!(eqs, D(v) ~ f / m + g)
    else
        push!(eqs, D(v) ~ f / m)
    end

    # position option
    if !isnothing(s_0)
        @parameters s_0 = s_0
        push!(pars, s_0)

        @variables s(t) = s_0
        push!(vars, s)

        push!(eqs, D(s) ~ v)
    end

    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [flange.v => v_0]),
                   flange)
end

const REL = Val(:relative)

"""
    Spring(; name, k, delta_s_0 = 0.0,  v_a_0=0.0, v_b_0=0.0)

Linear 1D translational spring

# Parameters:
- `k`: [N/m] Spring constant
- `delta_s_0`: initial spring stretch
- `v_a_0`: [m/s] Initial value of absolute linear velocity at flange_a (default 0 m/s)
- `v_b_0`: [m/s] Initial value of absolute linear velocity at flange_b (default 0 m/s)

# Connectors:
- `flange_a: 1-dim. translational flange on one side of spring`
- `flange_b: 1-dim. translational flange on opposite side of spring`
"""
function Spring(; name, k, delta_s_0 = 0.0, v_a_0 = 0.0, v_b_0 = 0.0)
    Spring(REL; name, k, delta_s_0, v_a_0, v_b_0)
end # default

function Spring(::Val{:relative}; name, k, delta_s_0 = 0.0, v_a_0 = 0.0, v_b_0 = 0.0)
    pars = @parameters begin
        k = k
        delta_s_0 = delta_s_0
        v_a_0 = v_a_0
        v_b_0 = v_b_0
    end
    vars = @variables begin
        delta_s(t) = delta_s_0
        f(t) = 0
    end

    @named flange_a = MechanicalPort()
    @named flange_b = MechanicalPort()

    eqs = [D(delta_s) ~ flange_a.v - flange_b.v
           f ~ k * delta_s
           flange_a.f ~ +f
           flange_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [flange_a.v => v_a_0, flange_b.v => v_b_0]), flange_a,
                   flange_b) #flange_a.f => +k*delta_s_0, flange_b.f => -k*delta_s_0
end

const ABS = Val(:absolute)
function Spring(::Val{:absolute}; name, k, s_a_0 = 0, s_b_0 = 0, v_a_0 = 0.0, v_b_0 = 0.0,
                l = 0)
    pars = @parameters begin
        k = k
        s_a_0 = s_a_0
        s_b_0 = s_b_0
        v_a_0 = v_a_0
        v_b_0 = v_b_0
        l = l
    end
    vars = @variables begin
        s1(t) = s_a_0
        s2(t) = s_b_0
        f(t) = 0
    end

    @named flange_a = MechanicalPort()
    @named flange_b = MechanicalPort()

    eqs = [D(s1) ~ flange_a.v
           D(s2) ~ flange_b.v
           f ~ k * (s1 - s2 - l) #delta_s
           flange_a.f ~ +f
           flange_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [flange_a.v => v_a_0, flange_b.v => v_b_0]), flange_a,
                   flange_b) #, flange_a.f => k * (s_a_0 - s_b_0 - l)
end

"""
    Damper(; name, d, v_a_0=0.0, v_b_0=0.0)

Linear 1D translational damper

# Parameters:
- `d`: [N.s/m] Damping constant
- `v_a_0`: [m/s] Initial value of absolute linear velocity at flange_a (default 0 m/s)
- `v_b_0`: [m/s] Initial value of absolute linear velocity at flange_b (default 0 m/s)

# Connectors:
- `flange_a: 1-dim. translational flange on one side of damper`
- `flange_b: 1-dim. translational flange on opposite side of damper`
"""
function Damper(; name, d, v_a_0 = 0.0, v_b_0 = 0.0)
    pars = @parameters begin
        d = d
        v_a_0 = v_a_0
        v_b_0 = v_b_0
    end
    vars = @variables begin
        v(t) = v_a_0 - v_b_0
        f(t) = 0.0
    end

    @named flange_a = MechanicalPort()
    @named flange_b = MechanicalPort()

    eqs = [v ~ flange_a.v - flange_b.v
           f ~ v * d
           flange_a.f ~ +f
           flange_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [flange_a.v => v_a_0, flange_b.v => v_b_0]), flange_a,
                   flange_b) #flange_a.f => +(v_a_0 - v_b_0)*d, flange_b.f => -(v_a_0 - v_b_0)*d
end
