@connector function Pin(;name)
    sts = @variables begin
        v(t)                    # Potential at the pin [V]
        i(t), [connect=Flow]    # Current flowing into the pin [A]
    end 
    ODESystem(Equation[], t, sts, [], name=name, defaults=Dict(v=>1.0, i=>1.0))
end
@doc """
    Pin(; name)

A pin in an analog circuit.

# States
- `v(t)`: [`V`] The voltage at this pin
- `i(t)`: [`A`] The current passing through this pin
""" Pin

"""
    OnePort(; name, v_start=0.0, i_start=0.0)

Component with two electrical pins `p` and `n` and current `i` from `p` to `n`.

# States
- `v(t)`: [`V`] The voltage across component `p.v - n.v`
- `i(t)`: [`A`] The current passing through positive pin

# Parameters:
- `v_start`: [`V`] Initial voltage across the component
- `i_start`: [`A`] Initial current through the component
"""
function OnePort(;name, v_start=0.0, i_start=0.0)
    @named p = Pin()
    @named n = Pin()
    sts = @variables begin
        v(t)=v_start
        i(t)=i_start
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
@doc """
    DigitalPin(; name)

A pin in a digital circuit.

# States
- `v(t)`: [`V`] The voltage at this pin
- `i(t)`: [`A`] The current passing through this pin
- `val(t)`: The binary value of the pin at this point. A voltage from `0V` to `0.8V` is a binary value of `0`. 
A voltage in the range `2.0V` to `5.0V` is `1`. Any other value is `X`.
""" DigitalPin

