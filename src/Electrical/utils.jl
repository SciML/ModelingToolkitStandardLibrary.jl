@connector function Pin(;name)
    @variables v(t)=0.0 i(t)=0.0 [connect = Flow]
    ODESystem(Equation[], t, [v, i], [], name=name)
end

@connector function DigitalPin(; name)
    @variables val(t)=0 v(t)=0.0 i(t)=0.0 [connect = Flow]
    eqs = [
        val ~ IfElse.ifelse((0.0 <= v) & (v <= 0.8) | (2.0 <= v) & (v <= 5.0),
                                IfElse.ifelse(v > 2.0, 1, 0), X)
    ]
    ODESystem(Equation[], t, [val, v, i], [], name=name)
end
