@connector function Pin(;name)
    @variables v(t) i(t)
    ODESystem(Equation[], t, [v, i], [], name=name, defaults=Dict(v=>1.0, i=>1.0))
end

@connector function DigitalPin(; name)
    @variables val(t) v(t) i(t)
    eqs = [
        val ~ IfElse.ifelse((0.0 <= v) & (v <= 0.8) | (2.0 <= v) & (v <= 5.0),
                                IfElse.ifelse(v > 2.0, 1, 0), X)
    ]
    ODESystem(Equation[], t, [val, v, i], [], defaults=Dict(val=>0, i=>0), name=name)
end

abstract type ElectricalPin end
ModelingToolkit.promote_connect_rule(::Type{DigitalPin}, ::Type{Pin}) = ElectricalPin
ModelingToolkit.promote_connect_rule(::Type{Pin}, ::Type{DigitalPin}) = ElectricalPin
ModelingToolkit.promote_connect_rule(::Type{ElectricalPin}, ::Type{DigitalPin}) = ElectricalPin
ModelingToolkit.promote_connect_rule(::Type{ElectricalPin}, ::Type{Pin}) = ElectricalPin

function ModelingToolkit.connect(::Type{<:Pin}, ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end

    return eqs
end

function ModelingToolkit.connect(::Type{DigitalPin}, ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].val ~ ps[i+1].val)
    end
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end
    return eqs
end

function ModelingToolkit.connect(::Type{ElectricalPin}, ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]

    # KVL
    digpins = ModelingToolkit.ODESystem[]
    for p in ps
        ModelingToolkit.get_connection_type(p) == DigitalPin && push!(digpins, p)
    end
    for i in 1:length(digpins)-1
        push!(eqs, digpins[i].val ~ digpins[i+1].val)
    end
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end
    return eqs
end
