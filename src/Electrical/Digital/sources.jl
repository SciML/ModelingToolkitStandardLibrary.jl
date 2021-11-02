function PulseDiff(; name, Val=1, dt=0.1)
    @named d = DigitalPin()
    @variables val(t)
    D = ModelingToolkit.Difference(t; dt=dt)
    
    eqs = [
        D(val) ~ Val
        val ~ d.val
    ]

    ODESystem(eqs, t, [val], [], systems=[d], defaults=Dict(Val=>0), name=name)
end

function Set(; name)
    @named d = DigitalPin()

    eqs = [
        d.val ~ 1
    ]
    ODESystem(eqs, t, [], [], systems=[d],  name=name)
end

function Reset(; name)
    @named d = DigitalPin()

    eqs = [
        d.val ~ 0
    ]
    ODESystem(eqs, t, [], [], systems=[d], name=name)
end

function Pulse(; name, duty_cycle=0.5, T=1.0)
    @named d = DigitalPin()

    eqs = [
        d.val ~ IfElse.ifelse(t%T > duty_cycle*T, 1, 0)
    ]
    ODESystem(eqs, t, [], [], systems=[d], name=name)
end
