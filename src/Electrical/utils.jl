@connector function Pin(;name)
    @variables v(t) i(t)
    ODESystem(Equation[], t, [v, i], [], name = name, defaults = Dict(v => 1.0, i => 1.0))
end

@doc """
```julia
@connector function Pin(; name)
```

A pin in an analog circuit.

# Variables
- `v(t)`
  The voltage at this pin
- `i(t)`
  The current passing through this pin
""" Pin

@connector function DigitalPin(; name)
    @variables val(t) v(t) i(t)
    eqs = [
        val ~ IfElse.ifelse((0.0 <= v) & (v <= 0.8) | (2.0 <= v) & (v <= 5.0),
                                IfElse.ifelse(v > 2.0, 1, 0), X)
    ]
    ODESystem(Equation[], t, [val, v, i], [], defaults = Dict(val => 0, i => 0), name = name)
end

@doc """
```julia
@connector function DigitalPin(; name)
```

A pin in a digital circuit.

# Variables
- `v(t)`
  The voltage at this pin
- `i(t)`
  The current passing through this pin
- `val(t)`
  The binary value of the pin at this point. A voltage from 0V to 0.8V is a binary value
  of 0. A voltage in the range 2.0V to 5.0V is 1. Any other value is X.
""" DigitalPin

abstract type ElectricalPin end
ModelingToolkit.promote_connect_rule(::Type{DigitalPin}, ::Type{Pin}) = ElectricalPin
ModelingToolkit.promote_connect_rule(::Type{Pin}, ::Type{DigitalPin}) = ElectricalPin
ModelingToolkit.promote_connect_rule(::Type{ElectricalPin}, ::Type{DigitalPin}) = ElectricalPin
ModelingToolkit.promote_connect_rule(::Type{ElectricalPin}, ::Type{Pin}) = ElectricalPin

"""
```julia
function ModelingToolkit.connect(::Type{<:Pin}, ps...)
function ModelingToolkit.connect(::Type{DigitalPin}, ps...)
function ModelingToolkit.connect(::Type{ElectricalPin}, ps...)
```

Returns equations for connecting the pins in `ps` of the specified type. Voltages
(and, in the case of `DigitalPin`s, values) of all provided pins are made equal, and
the total current flowing through them is 0.
"""
function ModelingToolkit.connect(::Type{<:Pin}, ps...)
    eqs = [
           0 ~ sum(p -> p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps) - 1
        push!(eqs, ps[i].v ~ ps[i + 1].v)
    end

    return eqs
end

function ModelingToolkit.connect(::Type{DigitalPin}, ps...)
    eqs = [
           0 ~ sum(p -> p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps) - 1
        push!(eqs, ps[i].val ~ ps[i + 1].val)
    end
    for i in 1:length(ps) - 1
        push!(eqs, ps[i].v ~ ps[i + 1].v)
    end
    return eqs
end

function ModelingToolkit.connect(::Type{ElectricalPin}, ps...)
    eqs = [
           0 ~ sum(p -> p.i, ps) # KCL
          ]

    # KVL
    digpins = ModelingToolkit.ODESystem[]
    for p in ps
        ModelingToolkit.get_connection_type(p) == DigitalPin && push!(digpins, p)
    end
    for i in 1:length(digpins) - 1
        push!(eqs, digpins[i].val ~ digpins[i + 1].val)
    end
    for i in 1:length(ps) - 1
        push!(eqs, ps[i].v ~ ps[i + 1].v)
    end
    return eqs
end
