"""
    Fixed(;name)

Fixes a port position (velocity = 0)

# Connectors:
- `port: 1-dim. translational port`
"""
function Fixed(; name)
    @named port = Port()
    eqs = [port.v ~ 0]
    return compose(ODESystem(eqs, t, [], []; name = name, defaults = [port.v => 0]), port)
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
- `port: 1-dim. translational port`
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

    @named port = Port()

    eqs = [port.v ~ v
           port.f ~ f]

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

    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [port.v => v_0]), port)
end

const REL = Val(:relative)

"""
    Spring(; name, k, delta_s_0 = 0.0,  v1_0=0.0, v2_0=0.0)

Linear 1D translational spring

# Parameters:
- `k`: [N/m] Spring constant
- `delta_s_0`: initial spring stretch
- `v1_0`: [m/s] Initial value of absolute linear velocity at port_a (default 0 m/s)
- `v2_0`: [m/s] Initial value of absolute linear velocity at port_b (default 0 m/s)

# Connectors:
- `port_a: 1-dim. translational port on one side of spring`
- `port_b: 1-dim. translational port on opposite side of spring`
"""
Spring(; name, k, delta_s_0 = 0.0, v1_0 = 0.0, v2_0 = 0.0) = Spring(REL; name, k, delta_s_0, v1_0, v2_0) # default 

function Spring(::Val{:relative}; name, k, delta_s_0 = 0.0, v1_0 = 0.0, v2_0 = 0.0)
    pars = @parameters begin
        k = k
        delta_s_0 = delta_s_0
        v1_0 = v1_0
        v2_0 = v2_0
    end
    vars = @variables begin
        delta_s(t) = delta_s_0
        f(t) = 0
    end

    @named port_a = Port()
    @named port_b = Port()

    eqs = [D(delta_s) ~ port_a.v - port_b.v
           f ~ k * delta_s
           port_a.f ~ +f
           port_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [port_a.v => v1_0, port_b.v => v2_0]), port_a, port_b) #port_a.f => +k*delta_s_0, port_b.f => -k*delta_s_0
end

const ABS = Val(:absolute)
function Spring(::Val{:absolute}; name, k, s1_0 = 0, s2_0 = 0, v1_0 = 0.0, v2_0 = 0.0, l = 0)
    pars = @parameters begin
        k = k
        s1_0 = s1_0
        s2_0 = s2_0
        v1_0 = v1_0
        v2_0 = v2_0
        l = l
    end
    vars = @variables begin
        s1(t) = s1_0
        s2(t) = s2_0
        f(t) = 0
    end

    @named port_a = Port()
    @named port_b = Port()

    eqs = [D(s1) ~ port_a.v
           D(s2) ~ port_b.v
           f ~ k * (s1 - s2 - l) #delta_s
           port_a.f ~ +f
           port_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [port_a.v => v1_0, port_b.v => v2_0]), port_a, port_b) #, port_a.f => k * (s1_0 - s2_0 - l)
end

"""
    Damper(; name, d, v1_0=0.0, v2_0=0.0)

Linear 1D translational damper

# Parameters:
- `d`: [N.s/m] Damping constant
- `v1_0`: [m/s] Initial value of absolute linear velocity at port_a (default 0 m/s)
- `v2_0`: [m/s] Initial value of absolute linear velocity at port_b (default 0 m/s)

# Connectors:
- `port_a: 1-dim. translational port on one side of damper`
- `port_b: 1-dim. translational port on opposite side of damper`
"""
function Damper(; name, d, v1_0 = 0.0, v2_0 = 0.0)
    pars = @parameters begin
        d = d
        v1_0 = v1_0
        v2_0 = v2_0
    end
    vars = @variables begin
        v(t) = v1_0 - v2_0
        f(t) = 0.0
    end

    @named port_a = Port()
    @named port_b = Port()

    eqs = [v ~ port_a.v - port_b.v
           f ~ v * d
           port_a.f ~ +f
           port_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [port_a.v => v1_0, port_b.v => v2_0]), port_a, port_b) #port_a.f => +(v1_0 - v2_0)*d, port_b.f => -(v1_0 - v2_0)*d
end
