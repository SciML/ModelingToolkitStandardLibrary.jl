@connector function HeatPort(; name)
    @variables T(t), Q_flow(t) # Temperature and Heat-flow-rate
    ODESystem(Equation[], t, [T, Q_flow], [], name=name)
end

function ModelingToolkit.connect(::Type{<:HeatPort}, ps...)
    eqs = [
           0 ~ sum(p->p.Q_flow, ps)
          ]
    for i in 1:length(ps)-1
        push!(eqs, ps[i].T ~ ps[i+1].T)
    end

    return eqs
end
