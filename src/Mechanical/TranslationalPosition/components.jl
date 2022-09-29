"""
    Fixed(;name, s₀=0.0)

Flange fixed in housing at a given position.

# Parameters:
- `s₀`: [m] Fixed offset position of housing

# Connectors:
- `T: 1-dim. translational flange`
"""
function Fixed(; name, s₀ = 0.0)
    pars = @parameters s₀ = s₀
    vars = []

    @named T = Flange()

    eqs = [T.s ~ s₀]

    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults=[T.s => s₀]), T)
end

"""
    Mass(; name, m, s₀ = 0.0, v₀ = 0.0)

Sliding mass with inertia

# Parameters:
- `m`: [kg] Mass of sliding mass
- `s₀`: [m] Initial value of absolute position of sliding mass
- `v₀`: [m/s] Initial value of absolute linear velocity of sliding mass

# States: 
- `s`: [m] Absolute position of sliding mass
- `v`: [m/s] Absolute linear velocity of sliding mass (= der(s)) 

# Connectors:
- `T: 1-dim. translational flange of mass`
"""
function Mass(; name, m, s₀ = 0.0, v₀ = 0.0)
    @named T = Flange()
    pars = @parameters begin
        m = m
        s₀=s₀
        v₀=v₀
    end 
    vars = @variables begin
        s(t) = s₀
        v(t) = v₀
        f(t) = 0
    end
    eqs = [T.s ~ s
           T.f ~ f
           D(s) ~ v
           D(v) ~ f / m]
    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults=[T.s => s₀]), T)
end


const REL = Val(:relative)
function Spring(::Val{:relative}; name, k, v1₀=0.0, v2₀=0.0, Δs0=0, s1₀ = 0, s2₀ = 0)
    pars = @parameters begin
        k=k 
        v1₀=v1₀
        v2₀=v2₀ 
        s1₀=s1₀ 
        s2₀=s2₀
        Δs0=Δs0
    end 
    vars = @variables begin
        v1(t) = v1₀
        v2(t) = v2₀ 
        Δs(t) = Δs0
        f(t) = Δs0*k
    end

    @named T1 = Flange()
    @named T2 = Flange()

    eqs = [D(T1.s) ~ v1
           D(T2.s) ~ v2
           D(Δs) ~ v1 - v2
           f ~ k * Δs
           T1.f ~ +f
           T2.f ~ -f]
    
    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [T1.s => s1₀, T2.s => s2₀, T1.f => +Δs0*k, T2.f => -Δs0*k]), T1, T2)
end

const ABS = Val(:absolute)

"""
    Spring(; name, k, s1₀ = 0, s2₀ = 0, l=0)

Linear 1D translational spring

# Parameters:
- `k`: [N/m] Spring constant
- `l`: Unstretched spring length
- `s1₀`: [m] Initial value of absolute position of T1
- `s2₀`: [m] Initial value of absolute position of T2

# Connectors:
- `T1: 1-dim. translational flange on one side of spring`
- `T2: 1-dim. translational flange on opposite side of spring`
"""
Spring(; name, k, s1₀ = 0, s2₀ = 0, l=0) = Spring(ABS; name, k, s1₀, s2₀, l) #default function

function Spring(::Val{:absolute}; name, k, s1₀ = 0, s2₀ = 0, l=0)
    pars = @parameters begin
        k=k 
        s1₀=s1₀ 
        s2₀=s2₀
        l=l
    end 
    vars = @variables begin
    # Δs(t) = s1 - s2
        f(t) = k * (s1₀ - s2₀ - l) 
    end

    @named T1 = Flange()
    @named T2 = Flange()

    eqs = [
           #    Δs ~ T1.s - T2.s
           f ~ k * (T1.s - T2.s - l) #Δs
           T1.f ~ +f
           T2.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [T1.s => s1₀, T2.s => s2₀, T1.f => k * (s1₀ - s2₀ - l), T2.f => -k * (s1₀ - s2₀ - l)]), T1, T2)
end

"""
    Damper(; name, d, v1₀=0.0, v2₀=0.0, s1₀ = 0, s2₀ = 0)

Linear 1D translational damper

# Parameters:
- `d`: [N.s/m] Damping constant
- `s1₀`: [m] Initial value of absolute position of T1
- `s2₀`: [m] Initial value of absolute position of T2
- `v1₀`: [m/s] Initial value of absolute linear velocity of T1
- `v2₀`: [m/s] Initial value of absolute linear velocity of T1

# Connectors:
- `T1: 1-dim. translational flange on one side of damper`
- `T2: 1-dim. translational flange on opposite side of damper`
"""
function Damper(; name, d, v1₀=0.0, v2₀=0.0, s1₀ = 0, s2₀ = 0)
    pars = @parameters begin
        d = d
        s1₀=s1₀ 
        s2₀=s2₀
        v1₀=v1₀
        v2₀=v2₀ 
    end
    vars = @variables begin
        v1(t) = v1₀
        v2(t) = v2₀
        f(t) = +(v1₀ - v2₀) * d
    end

    @named T1 = Flange()
    @named T2 = Flange()

    eqs = [D(T1.s) ~ v1
           D(T2.s) ~ v2
           f ~ (v1 - v2) * d
           T1.f ~ +f
           T2.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [T1.s => s1₀, T2.s => s2₀, T1.f => +(v1₀ - v2₀) * d, T2.f => -(v1₀ - v2₀) * d]), T1, T2)
end
