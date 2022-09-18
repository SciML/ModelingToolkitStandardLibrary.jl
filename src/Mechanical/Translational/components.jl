function Fixed(; name)
    @named port = Port()
    eqs = [port.v ~ 0]
    return compose(ODESystem(eqs, t, [], []; name = name), port)
end

function Body(; name, v0 = 0.0, m = 1.0, s0=nothing, g=nothing)
    @named port = Port(v0 = v0)
    pars = @parameters m=m
    vars = @variables begin
        v(t) = v0
        f(t) = m*v0
    end

    eqs = [
        port.v ~ v
        port.f ~ f
    ]

    # gravity option
    if !isnothing(g)
        @parameters g=g
        push!(pars, g)
        push!(eqs, D(v) ~ f/m+g)
    else
        push!(eqs, D(v) ~ f/m)
    end

    # position option
    if !isnothing(s0)
        @parameters s0=s0
        push!(pars, s0)

        @variables s(t)=s0
        push!(vars, s)
        
        push!(eqs, D(s) ~ v)
    end
    
    return compose(ODESystem(eqs, t, vars, pars; name = name), port)
end


function Spring(;name, k = 1e3, delta0 = 0.0)
    pars = @parameters k=k
    vars = @variables begin
        Δx(t) = delta0
        f(t) = k*delta0
    end
    
    @named port_a = Port()
    @named port_b = Port()
    
    eqs = [
        D(Δx) ~ port_a.v - port_b.v
        f ~ k*Δx
        port_a.f ~ +f
        port_b.f ~ -f
    ]
    return compose(ODESystem(eqs, t, vars, pars; name = name), port_a, port_b)
end

function Damper(; name, d=1e2)
    pars = @parameters d=d
    vars = @variables v(t)=0.0 f(t)=0.0
    
    @named port_a = Port()
    @named port_b = Port()
    
    eqs = [
        v ~ port_a.v - port_b.v
        f ~ v*d

        port_a.f ~ +f
        port_b.f ~ -f
    ]
    return compose(ODESystem(eqs, t, vars, pars; name = name), port_a, port_b)
end


