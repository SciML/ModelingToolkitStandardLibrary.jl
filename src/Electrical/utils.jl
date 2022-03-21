@connector function Pin(;name)
    sts = @variables begin
        v(t)                    # Potential at the pin [V]
        i(t), [connect=Flow]    # Current flowing into the pin [A]
    end 
    ODESystem(Equation[], t, sts, [], name=name, defaults=Dict(v=>1.0, i=>1.0))
end

function OnePort(;name, 
    v0=0.0, # [V] Initial voltage across the component
    i0=0.0, # [A] Initial current through the component
    )

    @named p = Pin()
    @named n = Pin()
    sts = @variables begin
        v(t)=v0
        i(t)=i0
    end
    eqs = [
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
    ]
    
    return compose(ODESystem(eqs, t, sts, []; name=name), p, n)
end

@connector function DigitalPin(; name)
    @variables val(t) v(t) i(t)
    eqs = [
        val ~ IfElse.ifelse((0.0 <= v) & (v <= 0.8) | (2.0 <= v) & (v <= 5.0),
                                IfElse.ifelse(v > 2.0, 1, 0), X)
    ]
    ODESystem(Equation[], t, [val, v, i], [], defaults=Dict(val=>0, i=>0), name=name)
end

