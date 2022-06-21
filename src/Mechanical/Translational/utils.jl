@connector function Flange(; name)
    sts = @variables begin
        s(t)
        f(t), [connect = Flow]
    end
    ODESystem(Equation[], t, sts, [], name=name, defaults=Dict(s => 0.0, f => 0.0))
end

@connector function Support(; name)
    sts = @variables begin
        s(t)
        f(t), [connect = Flow]
    end
    ODESystem(Equation[], t, sts, [], name=name, defaults=Dict(s => 0.0, f => 0.0))
end

function PartialCompliant(; name, s_rel_start=0.0, f_start=0.0)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sts = @variables begin
        s_rel(t) = s_rel_start
        f(t) = f_start
    end
    eqs = [
        s_rel ~ flange_b.s - flange_a.s
        flange_b.f ~ f
        flange_a.f ~ -f
    ]
    return compose(ODESystem(eqs, t, sts, []; name=name), flange_a, flange_b)
end

function PartialCompliantWithRelativeStates(; name, s_rel_start=0.0, v_start=0.0, a_start=0.0, f_start=0.0)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sts = @variables begin
        s_rel(t) = s_rel_start
        v_rel(t) = v_start
        a_rel(t) = a_start
        f(t) = f_start
    end
    eqs = [
        s_rel ~ flange_b.s - flange_a.s
        D(s_rel) ~ v_rel
        D(v_rel) ~ a_rel
        flange_b.f ~ f
        flange_a.f ~ -f
    ]
    return compose(ODESystem(eqs, t, sts, []; name=name), flange_a, flange_b)
end

function PartialElementaryOneFlangeAndSupport2(; name, use_support=false)
    @named flange = Flange()
    sys = [flange]
    @variables s_support(t)
    if use_support
        @named support = Support()
        eqs = [
            support.s ~ s_support
            support.f ~ -flange.f
        ]
        push!(sys, support)
    else
        eqs = [s_support ~ 0]
    end
    return compose(ODESystem(eqs, t, [s_support], []; name=name), sys)
end

function PartialElementaryTwoFlangesAndSupport2(; name, use_support=false)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sys = [flange_a, flange_b]
    @variables s_support(t) = 0.0
    if use_support
        @named support = Support()
        eqs = [
            support.s ~ s_support
            support.f ~ -flange_a.f - flange_b.f
        ]
        push!(sys, support)
    else
        eqs = [s_support ~ 0]
    end
    return compose(ODESystem(eqs, t, [s_support], []; name=name), sys)
end
