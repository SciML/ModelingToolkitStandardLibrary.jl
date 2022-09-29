"""
    Fixed(;name)

Fixes a port position (velocity = 0)

# Connectors:
- `T: 1-dim. translational port`
"""
function Fixed(; name)
    @named T = Port()
    eqs = [T.v ~ 0]
    return compose(ODESystem(eqs, t, [], []; name = name, defaults = [T.v => 0]), T)
end

"""
    Mass(; name, v₀ = 0.0, m, s₀ = nothing, g = nothing)

Sliding mass with inertia

# Parameters:
- `m`: [kg] mass of sliding body
- `v₀`: [m/s] Initial value of absolute linear velocity of sliding mass (default 0 m/s)
- `s₀`: [m] (optional) initial value of absolute position of sliding mass
- `g`: [m/s²] (optional) gravity field acting on the mass, positive value acts in the positive direction


# States: 
- `v`: [m/s] absolute linear velocity of sliding mass
- `s`: [m] (optional with parameter s₀) absolute position of sliding mass


# Connectors:
- `T: 1-dim. translational port`
"""
function Mass(; name, v₀ = 0.0, m, s₀ = nothing, g = nothing)
    pars = @parameters begin
        m = m
        v₀ = v₀
    end
    vars = @variables begin
        v(t) = v₀
        f(t) = 0
    end

    @named T = Port()

    eqs = [T.v ~ v
           T.f ~ f]

    # gravity option
    if !isnothing(g)
        @parameters g = g
        push!(pars, g)
        push!(eqs, D(v) ~ f / m + g)
    else
        push!(eqs, D(v) ~ f / m)
    end

    # position option
    if !isnothing(s₀)
        @parameters s₀ = s₀
        push!(pars, s₀)

        @variables s(t) = s₀
        push!(vars, s)

        push!(eqs, D(s) ~ v)
    end

    return compose(ODESystem(eqs, t, vars, pars; name = name, defaults = [T.v => v₀]), T)
end

const REL = Val(:relative)

"""
    Spring(; name, k, Δs₀ = 0.0,  v1₀=0.0, v2₀=0.0)

Linear 1D translational spring

# Parameters:
- `k`: [N/m] Spring constant
- `Δs₀`: initial spring stretch
- `v1₀`: [m/s] Initial value of absolute linear velocity at T1 (default 0 m/s)
- `v2₀`: [m/s] Initial value of absolute linear velocity at T2 (default 0 m/s)

# Connectors:
- `T1: 1-dim. translational port on one side of spring`
- `T2: 1-dim. translational port on opposite side of spring`
"""
Spring(; name, k, Δs₀ = 0.0, v1₀ = 0.0, v2₀ = 0.0) = Spring(REL; name, k, Δs₀, v1₀, v2₀) # default 

function Spring(::Val{:relative}; name, k, Δs₀ = 0.0, v1₀ = 0.0, v2₀ = 0.0)
    pars = @parameters begin
        k = k
        Δs₀ = Δs₀
        v1₀ = v1₀
        v2₀ = v2₀
    end
    vars = @variables begin
        Δs(t) = Δs₀
        f(t) = 0
    end

    @named T1 = Port()
    @named T2 = Port()

    eqs = [D(Δs) ~ T1.v - T2.v
           f ~ k * Δs
           T1.f ~ +f
           T2.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [T1.v => v1₀, T2.v => v2₀]), T1, T2) #T1.f => +k*Δs₀, T2.f => -k*Δs₀
end

const ABS = Val(:absolute)
function Spring(::Val{:absolute}; name, k, s1₀ = 0, s2₀ = 0, v1₀ = 0.0, v2₀ = 0.0, l = 0)
    pars = @parameters begin
        k = k
        s1₀ = s1₀
        s2₀ = s2₀
        v1₀ = v1₀
        v2₀ = v2₀
        l = l
    end
    vars = @variables begin
        s1(t) = s1₀
        s2(t) = s2₀
        f(t) = 0
    end

    @named T1 = Port()
    @named T2 = Port()

    eqs = [D(s1) ~ T1.v
           D(s2) ~ T2.v
           f ~ k * (s1 - s2 - l) #Δs
           T1.f ~ +f
           T2.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [T1.v => v1₀, T2.v => v2₀]), T1, T2) #, T1.f => k * (s1₀ - s2₀ - l)
end

"""
    Damper(; name, d, v1₀=0.0, v2₀=0.0)

Linear 1D translational damper

# Parameters:
- `d`: [N.s/m] Damping constant
- `v1₀`: [m/s] Initial value of absolute linear velocity at T1 (default 0 m/s)
- `v2₀`: [m/s] Initial value of absolute linear velocity at T2 (default 0 m/s)

# Connectors:
- `T1: 1-dim. translational port on one side of damper`
- `T2: 1-dim. translational port on opposite side of damper`
"""
function Damper(; name, d, v1₀ = 0.0, v2₀ = 0.0)
    pars = @parameters begin
        d = d
        v1₀ = v1₀
        v2₀ = v2₀
    end
    vars = @variables begin
        v(t) = v1₀ - v2₀
        f(t) = 0.0
    end

    @named T1 = Port()
    @named T2 = Port()

    eqs = [v ~ T1.v - T2.v
           f ~ v * d
           T1.f ~ +f
           T2.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name,
                             defaults = [T1.v => v1₀, T2.v => v2₀]), T1, T2) #T1.f => +(v1₀ - v2₀)*d, T2.f => -(v1₀ - v2₀)*d
end
