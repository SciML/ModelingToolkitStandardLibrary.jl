function Fixed(; name)
    @named T = Port()
    eqs = [T.v ~ 0]
    return compose(ODESystem(eqs, t, [], []; name = name), T)
end

function Mass(; name, v0 = 0.0, m, s0 = nothing, g = nothing)
    @named T = Port(v0 = v0)
    pars = @parameters m = m
    vars = @variables begin
        v(t) = v0
        f(t) = 0
    end

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
    if !isnothing(s0)
        @parameters s0 = s0
        push!(pars, s0)

        @variables s(t) = s0
        push!(vars, s)

        push!(eqs, D(s) ~ v)
    end

    return compose(ODESystem(eqs, t, vars, pars; name = name), T)
end

function Spring(; name, k, Δs0 = 0.0)
    pars = @parameters k=k Δs0=Δs0
    vars = @variables begin
        Δs(t) = Δs0
        f(t) = 0
    end

    @named T1 = Port()
    @named T2 = Port()

    eqs = [D(Δs) ~ T1.v - T2.v
           f ~ k * Δs
           T1.f ~ +f
           T2.f ~ -f]
    return ODESystem(eqs, t, vars, pars; name = name, systems = [T1, T2])
end

function Damper(; name, d)
    pars = @parameters d = d
    vars = @variables v(t)=0.0 f(t)=0.0

    @named T1 = Port()
    @named T2 = Port()

    eqs = [v ~ T1.v - T2.v
           f ~ v * d
           T1.f ~ +f
           T2.f ~ -f]
    return compose(ODESystem(eqs, t, vars, pars; name = name), T1, T2)
end
