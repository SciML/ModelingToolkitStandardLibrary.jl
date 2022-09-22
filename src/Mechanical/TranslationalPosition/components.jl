"""
    Fixed(;name, s0=0.0)

Flange fixed in housing at a given position.

# Parameters:
- `s0`: [m] Fixed offset position of housing

# Connectors:
- `flange: 1-dim. translational flange`
"""
function Fixed(; name, s0 = 0.0)
    @named T = Flange()
    @parameters s0 = s0
    eqs = [T.s ~ s0]
    return compose(ODESystem(eqs, t, [], [s0]; name = name), T)
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
function Mass(; name, m, s0 = 0.0, v0 = 0.0)
    @named T = Flange()
    pars = @parameters m = m
    vars = @variables begin
        s(t) = s0
        v(t) = v0
        f(t) = 0
    end
    eqs = [T.s ~ s
           T.f ~ f
           D(s) ~ v
           D(v) ~ f / m]
    return compose(ODESystem(eqs, t, vars, pars; name = name), T)
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
const REL = Val(:relative)
function Spring(::Val{:relative}; name, k, Δs0 = 0.0)
    pars = @parameters k=k Δs0=Δs0
    vars = @variables begin
        v1(t) = 0
        v2(t) = 0
        Δs(t) = Δs0
        f(t) = 0
    end

    @named T1 = Flange()
    @named T2 = Flange()

    eqs = [D(T1.s) ~ v1
           D(T2.s) ~ v2
           D(Δs) ~ v1 - v2
           f ~ k * Δs
           T1.f ~ +f
           T2.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name), T1, T2)
end

const ABS = Val(:absolute)
Spring(; name, k, s1₀ = 0, s2₀ = 0) = Spring(ABS; name, k, s1₀, s2₀) #default function
function Spring(::Val{:absolute}; name, k, s1₀ = 0, s2₀ = 0)
    pars = @parameters k=k s1=s1₀ s2=s2₀
    vars = @variables begin
    # Δs(t) = s1 - s2
    f(t) = 0 end

    @named T1 = Flange()
    @named T2 = Flange()

    eqs = [
           #    Δs ~ T1.s - T2.s
           f ~ k * (T1.s - T2.s) #Δs
           T1.f ~ +f
           T2.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [T1.s => s1, T2.s => s2]), T1, T2)
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
    pars = @parameters d = d
    vars = @variables begin
        v1(t) = 0
        v2(t) = 0
        f(t) = 0
    end

    @named T1 = Flange()
    @named T2 = Flange()

    eqs = [D(T1.s) ~ v1
           D(T2.s) ~ v2
           f ~ (v1 - v2) * d
           T1.f ~ +f
           T2.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [T1.s => 0, T2.s => 0]), T1, T2)
end
