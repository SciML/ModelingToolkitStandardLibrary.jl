
@enum BoundaryType NoBoundary SpringDamperBoundary ForcedBoundary

#Note: boundary is an argument because this is a compile time constant
function Body(boundary = NoBoundary; name, x0 = 0.0, v0 = 0.0, mass = 1.0,  g=0.0, upper_limit = +Inf, lower_limit = -Inf, k=1e9, d=1e6)
    @named T = Port(v0 = v0, f0 = mass*g)
    @parameters mass=mass g=g
    @variables begin
        x(t) = x0
        dx(t) = v0
        ddx(t) = 0.0
    end

    if boundary != NoBoundary
        @variables begin
            f_floor(t) = 0.0
            f_ceil(t) = 0.0
        end
        @parameters begin
            ubnd = upper_limit
            lbnd = lower_limit
            k=k
            d=d
        end
    end

    eqs = [
        T.v ~ dx

        D(x) ~ dx
        D(dx) ~ ddx
    ]

    if boundary == NoBoundary
        feqs = [
            ddx * mass ~ mass*g + T.f
        ]
        return compose(ODESystem([eqs; feqs], t, [x, dx, ddx], [mass, g]; name = name), T)
    elseif boundary == SpringDamperBoundary

        # let
        Δceil = x - ubnd
        Δfloor = lbnd - x

        feqs = [
            ddx * mass ~ mass*g + T.f - f_floor - f_ceil
            f_ceil ~ ifelse(x < ubnd, 0.0, Δceil*k + dx*d)
            f_floor ~ ifelse(x > lbnd, 0.0, Δfloor*k + dx*d)
        ]
        return compose(ODESystem([eqs; feqs], t, [x, dx, ddx, f_ceil, f_floor], [mass, g, ubnd, lbnd, k, d]; name = name), T)
    elseif boundary == ForcedBoundary

        # let 
        f = mass*g + T.f # external forces

        feqs = [
            ddx * mass ~ f - f_floor - f_ceil
            0 ~ ifelse((x >= ubnd), 
                    ( x ) - ( ubnd ), # make x=>ubnd, solve for f_ceil to make this true
                    ( f_ceil ) - ( 0 )
                )
            0 ~ ifelse((x <= lbnd) , 
                ( x ) - ( lbnd ), # make x=>lbnd, solve for f_floor to make this true
                ( f_floor ) - ( 0 )
            )
        ]
        return compose(ODESystem([eqs; feqs], t, [x, dx, ddx, f_ceil, f_floor], [mass, g, ubnd, lbnd]; name = name), T)
    end
        
    
    
end


function Spring(;name, k = 1e3, delta0 = 0.0)
    @parameters k=k
    @variables begin
        x(t) = delta0
        dx(t) = 0.0
        f(t) 
    end
    
    @named T1 = Port()
    @named T2 = Port()
    
    eqs = [
        D(x) ~ dx
        dx ~ T1.v - T2.v
        f ~ k*x
        T1.f ~ +f
        T2.f ~ -f
    ]
    return compose(ODESystem(eqs, t, [x, dx, f], [k]; name = name), T1, T2)
end

function Damper(; name, c=1e2)
    @parameters c=c
    @variables dx=0.0
    
    @named T1 = Port()
    @named T2 = Port()
    
    eqs = [
        T1.v - T2.v ~ dx
        T1.f + T2.f ~ dx*c
    ]
    return compose(ODESystem(eqs, t, [x, dx, f], [k]; name = name), T1, T2)
end




