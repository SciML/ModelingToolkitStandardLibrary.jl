@connector function Pin(; name)
    sts = @variables begin
        v(t)                    # Potential at the pin [V]
        i(t), [connect = Flow]    # Current flowing into the pin [A]
    end
    ODESystem(Equation[], t, sts, [], name = name, defaults = Dict(v => 1.0, i => 1.0))
end
@doc """
    Pin(; name)

A pin in an analog circuit.

# States:
- `v(t)`: [`V`] The voltage at this pin
- `i(t)`: [`A`] The current passing through this pin
""" Pin

"""
    OnePort(; name, v_start=0.0, i_start=0.0)

Component with two electrical pins `p` and `n` and current `i` flows from `p` to `n`.

# States:
- `v(t)`: [`V`] The voltage across component `p.v - n.v`
- `i(t)`: [`A`] The current passing through positive pin

# Parameters:
- `v_start`: [`V`] Initial voltage across the component
- `i_start`: [`A`] Initial current through the component

# Connectors:
- `p` Positive pin
- `n` Negative pin
"""
function OnePort(; name, v_start = 0.0, i_start = 0.0)
    @named p = Pin()
    @named n = Pin()
    sts = @variables begin
        v(t) = v_start
        i(t) = i_start
    end
    eqs = [v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i]
    return compose(ODESystem(eqs, t, sts, []; name = name), p, n)
end

"""
    TwoPort(; name, v1_start=0.0, v2_start=0.0, i1_start=0.0, i2_start=0.0)

Component with four electrical pins `p1`, `n1`, `p2` and `n2`
Current `i1` flows from `p1` to `n1` and `i2` from `p2` to `n2`.

# States:
- `v1(t)`: [`V`] The voltage across first ports `p1.v - n1.v`
- `v2(t)`: [`V`] The voltage across second ports `p2.v - n2.v`
- `i1(t)`: [`A`] The current passing through positive pin `p1`
- `i2(t)`: [`A`] The current passing through positive pin `p2`

# Parameters:
- `v1_start`: [`V`] Initial voltage across p1 and n1
- `v2_start`: [`V`] Initial voltage across p2 and n2
- `i2_start`: [`A`] Initial current through p1
- `i2_start`: [`A`] Initial current through p2

# Connectors:
- `p1` First positive pin
- `p2` Second positive pin
- `n1` First negative pin
- `n2` Second Negative pin
"""

function TwoPort(; name, v1_start = 0.0, v2_start = 0.0, i1_start = 0.0, i2_start = 0.0)
    @named p1 = Pin()
    @named n1 = Pin()
    @named p2 = Pin()
    @named n2 = Pin()
    sts = @variables begin
        v1(t) = v1_start
        i1(t) = i1_start
        v2(t) = v2_start
        i2(t) = i2_start
    end
    eqs = [v1 ~ p1.v - n1.v
           0 ~ p1.i + n1.i
           i1 ~ p1.i
           v2 ~ p2.v - n2.v
           0 ~ p2.i + n2.i
           i2 ~ p2.i]
    return compose(ODESystem(eqs, t, sts, []; name = name), p1, p2, n1, n2)
end

@connector function DigitalPin(; name)
    @variables val(t) v(t) i(t)
    eqs = [
        val ~ IfElse.ifelse((0.0 <= v) & (v <= 0.8) | (2.0 <= v) & (v <= 5.0),
                            IfElse.ifelse(v > 2.0, 1, 0), X),
    ]
    ODESystem(Equation[], t, [val, v, i], [], defaults = Dict(val => 0, i => 0),
              name = name)
end
@doc """
    DigitalPin(; name)

A pin in a digital circuit.

# States:
- `v(t)`: [`V`] The voltage at this pin
- `i(t)`: [`A`] The current passing through this pin
- `val(t)`: The binary value of the pin at this point. A voltage from `0V` to `0.8V` is a binary value of `0`.
A voltage in the range `2.0V` to `5.0V` is `1`. Any other value is `X`.
""" DigitalPin
