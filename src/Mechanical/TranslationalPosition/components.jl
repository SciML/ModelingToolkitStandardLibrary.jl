"""
    Fixed(;name, s_0=0.0)

Flange fixed in housing at a given position.

# Parameters:
- `s_0`: [m] Fixed offset position of housing

# Connectors:
- `port: 1-dim. translational flange`
"""
function Fixed(; name, s_0 = 0.0)
    pars = @parameters s_0 = s_0
    vars = []

    @named port = Flange()

    eqs = [port.s ~ s_0]

    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [port.s => s_0]),
                   port)
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
- `port: 1-dim. translational flange of mass`
"""
function Mass(; name, m, s_0 = 0.0, v_0 = 0.0)
    @named port = Flange()
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
    eqs = [port.s ~ s
           port.f ~ f
           D(s) ~ v
           D(v) ~ f / m]
    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [port.s => s_0]),
                   port)
end

const REL = Val(:relative)
function Spring(::Val{:relative}; name, k, v1_0 = 0.0, v2_0 = 0.0, delta_s0 = 0, s1_0 = 0,
                s2_0 = 0)
    pars = @parameters begin
        k = k
        v1_0 = v1_0
        v2_0 = v2_0
        s1_0 = s1_0
        s2_0 = s2_0
        delta_s0 = delta_s0
    end
    vars = @variables begin
        v1(t) = v1_0
        v2(t) = v2_0
        delta_s(t) = delta_s0
        f(t) = delta_s0 * k
    end

    @named port_a = Flange()
    @named port_b = Flange()

    eqs = [D(port_a.s) ~ v1
           D(port_b.s) ~ v2
           D(delta_s) ~ v1 - v2
           f ~ k * delta_s
           port_a.f ~ +f
           port_b.f ~ -f]

    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [
                                 port_a.s => s1_0,
                                 port_b.s => s2_0,
                                 port_a.f => +delta_s0 * k,
                                 port_b.f => -delta_s0 * k,
                             ]), port_a, port_b)
end

const ABS = Val(:absolute)

"""
    Spring(; name, k, s1_0 = 0, s2_0 = 0, l=0)

Linear 1D translational spring

# Parameters:
- `k`: [N/m] Spring constant
- `l`: Unstretched spring length
- `s1_0`: [m] Initial value of absolute position of port_a
- `s2_0`: [m] Initial value of absolute position of port_b

# Connectors:
- `port_a: 1-dim. translational flange on one side of spring`
- `port_b: 1-dim. translational flange on opposite side of spring`
"""
Spring(; name, k, s1_0 = 0, s2_0 = 0, l = 0) = Spring(ABS; name, k, s1_0, s2_0, l) #default function

function Spring(::Val{:absolute}; name, k, s1_0 = 0, s2_0 = 0, l = 0)
    pars = @parameters begin
        k = k
        s1_0 = s1_0
        s2_0 = s2_0
        l = l
    end
    vars = @variables begin
    # delta_s(t) = s1 - s2
    f(t) = k * (s1_0 - s2_0 - l) end

    @named port_a = Flange()
    @named port_b = Flange()

    eqs = [
           #    delta_s ~ port_a.s - port_b.s
           f ~ k * (port_a.s - port_b.s - l) #delta_s
           port_a.f ~ +f
           port_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [
                                 port_a.s => s1_0,
                                 port_b.s => s2_0,
                                 port_a.f => k * (s1_0 - s2_0 - l),
                                 port_b.f => -k * (s1_0 - s2_0 - l),
                             ]), port_a, port_b)
end

"""
    Damper(; name, d, v1_0=0.0, v2_0=0.0, s1_0 = 0, s2_0 = 0)

Linear 1D translational damper

# Parameters:
- `d`: [N.s/m] Damping constant
- `s1_0`: [m] Initial value of absolute position of port_a
- `s2_0`: [m] Initial value of absolute position of port_b
- `v1_0`: [m/s] Initial value of absolute linear velocity of port_a
- `v2_0`: [m/s] Initial value of absolute linear velocity of port_a

# Connectors:
- `port_a: 1-dim. translational flange on one side of damper`
- `port_b: 1-dim. translational flange on opposite side of damper`
"""
function Damper(; name, d, v1_0 = 0.0, v2_0 = 0.0, s1_0 = 0, s2_0 = 0)
    pars = @parameters begin
        d = d
        s1_0 = s1_0
        s2_0 = s2_0
        v1_0 = v1_0
        v2_0 = v2_0
    end
    vars = @variables begin
        v1(t) = v1_0
        v2(t) = v2_0
        f(t) = +(v1_0 - v2_0) * d
    end

    @named port_a = Flange()
    @named port_b = Flange()

    eqs = [D(port_a.s) ~ v1
           D(port_b.s) ~ v2
           f ~ (v1 - v2) * d
           port_a.f ~ +f
           port_b.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [
                                 port_a.s => s1_0,
                                 port_b.s => s2_0,
                                 port_a.f => +(v1_0 - v2_0) * d,
                                 port_b.f => -(v1_0 - v2_0) * d,
                             ]), port_a, port_b)
end
