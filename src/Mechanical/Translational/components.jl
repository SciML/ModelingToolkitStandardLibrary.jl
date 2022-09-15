function Fixed(; name)
    @named port = Port()
    eqs = [port.v ~ 0]
    return compose(ODESystem(eqs, t, [], []; name = name), port)
end

function Body(; name, x0 = 0.0, v0 = 0.0, m = 1.0,  g=0.0)
    @named port = Port(v0 = v0, f0 = m*g)
    @parameters m=m g=g
    @variables begin
        x(t) = x0
        dx(t) = v0
        ddx(t) = 0.0
    end

    eqs = [
        port.v ~ dx

        D(x) ~ dx
        D(dx) ~ ddx

        ddx * m ~ m*g + port.f
    ]
    
    return compose(ODESystem(eqs, t, [x, dx, ddx], [m, g]; name = name), port)
end


function Spring(;name, k = 1e3, delta0 = 0.0)
    @parameters k=k
    @variables begin
        x(t) = delta0
        dx(t) = 0.0
        f(t) = k*x
    end
    
    @named port_a = Port()
    @named port_b = Port()
    
    eqs = [
        D(x) ~ dx
        dx ~ port_a.v - port_b.v
        f ~ k*x
        port_a.f ~ +f
        port_b.f ~ -f
    ]
    return compose(ODESystem(eqs, t, [x, dx, f], [k]; name = name), port_a, port_b)
end

function Damper(; name, d=1e2)
    @parameters d=d
    @variables dx(t)=0.0 f(t)=0.0
    
    @named port_a = Port()
    @named port_b = Port()
    
    eqs = [
        port_a.v - port_b.v ~ dx
        f ~ dx*d
        port_a.f ~ +f
        port_b.f ~ -f
    ]
    return compose(ODESystem(eqs, t, [dx, f], [d]; name = name), port_a, port_b)
end


