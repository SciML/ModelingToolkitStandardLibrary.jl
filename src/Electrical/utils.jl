@connector function Pin(; name, v = nothing, i = nothing)
    vars = @variables begin
        v(t) = v                  # Potential at the pin [V]
        i(t) = i, [connect = Flow]    # Current flowing into the pin [A]
    end
    System(Equation[], t, vars, []; name)
end
@doc """
    Pin(; name)

A pin in an analog circuit.

# States:
- `v(t)`: [`V`] The voltage at this pin
- `i(t)`: [`A`] The current passing through this pin
""" Pin

"""
    OnePort(; name, v = nothing, i = nothing)

Component with two electrical pins `p` and `n` and current `i` flows from `p` to `n`.

# States:

  - `v(t)`: [`V`] The voltage across component `p.v - n.v`
  - `i(t)`: [`A`] The current passing through positive pin

# Connectors:

  - `p` Positive pin
  - `n` Negative pin
"""
@component function OnePort(; v = nothing, i = nothing, name)
    pars = @parameters begin
    end

    systems = @named begin
        p = Pin()
        n = Pin()
    end

    vars = @variables begin
        v(t) = v
        i(t) = i
    end

    equations = Equation[
        v ~ p.v - n.v,
        0 ~ p.i + n.i,
        i ~ p.i
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    TwoPort(; name, v1 = nothing, v2 = nothing, i1 = nothing, i2 = nothing)

Component with four electrical pins `p1`, `n1`, `p2` and `n2`
Current `i1` flows from `p1` to `n1` and `i2` from `p2` to `n2`.

# States:
- `v1(t)`: [`V`] The voltage across first ports `p1.v - n1.v`
- `v2(t)`: [`V`] The voltage across second ports `p2.v - n2.v`
- `i1(t)`: [`A`] The current passing through positive pin `p1`
- `i2(t)`: [`A`] The current passing through positive pin `p2`

# Connectors:
- `p1` First positive pin
- `p2` Second positive pin
- `n1` First negative pin
- `n2` Second Negative pin
"""

@component function TwoPort(; v1 = nothing, v2 = nothing, i1 = nothing, i2 = nothing, name)
    pars = @parameters begin
    end

    systems = @named begin
        p1 = Pin()
        n1 = Pin()
        p2 = Pin()
        n2 = Pin()
    end

    vars = @variables begin
        v1(t) = v1
        i1(t) = i1
        v2(t) = v2
        i2(t) = i2
    end

    equations = Equation[
        v1 ~ p1.v - n1.v,
        0 ~ p1.i + n1.i,
        i1 ~ p1.i,
        v2 ~ p2.v - n2.v,
        0 ~ p2.i + n2.i,
        i2 ~ p2.i
    ]

    return System(equations, t, vars, pars; name, systems)
end

@connector function DigitalPin(; name)
    @variables val(t) v(t) i(t)
    eqs = [
        val ~ IfElse.ifelse((0.0 <= v) & (v <= 0.8) | (2.0 <= v) & (v <= 5.0),
        IfElse.ifelse(v > 2.0, 1, 0), X)
    ]
    System(Equation[], t, [val, v, i], [], guesses = Dict(val => 0, i => 0),
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
