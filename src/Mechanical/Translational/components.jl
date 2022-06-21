function Fixed(; name, s0=0.0)
    @named flange = Flange()
    @parameters s0 = s0
    eqs = [flange.s ~ s0]
    return compose(ODESystem(eqs, t, [], [s0]; name=name), flange)
end

function Mass(; name, m, s_start=0.0, v_start=0.0, a_start=0.0)
    @named flange_a = Flange()
    @named flange_b = Flange()
    @parameters m = m
    sts = @variables begin
        s(t) = s_start
        v(t) = v_start
        a(t) = a_start
    end
    eqs = [
        s ~ flange_a.s
        s ~ flange_b.s
        D(s) ~ v
        D(v) ~ a
        m * a ~ flange_a.f + flange_b.f
    ]
    return compose(ODESystem(eqs, t, sts, [m]; name=name), flange_a, flange_b)
end

function Spring(; name, c, s_rel0=0.0)
    @named partial_comp = PartialCompliant()
    @unpack s_rel, f = partial_comp
    pars = @parameters begin
        c = c
        s_rel0 = s_rel0
    end
    eqs = [f ~ c * (s_rel - s_rel0)]
    extend(ODESystem(eqs, t, [], pars; name=name), partial_comp)
end

function Damper(; name, d)
    @named partial_comp = PartialCompliantWithRelativeStates()
    @unpack v_rel, f = partial_comp
    pars = @parameters d = d
    eqs = [f ~ d * v_rel]
    extend(ODESystem(eqs, t, [], pars; name=name), partial_comp)
end

function IdealGear(; name, ratio, use_support=false)
    @named partial_element = PartialElementaryTwoFlangesAndSupport2(use_support=use_support)
    @unpack s_support, flange_a, flange_b = partial_element
    @parameters ratio = ratio
    sts = @variables s_a(t) = 0.0 s_b(t) = 0.0
    eqs = [
        s_a ~ flange_a.s - s_support
        s_b ~ flange_b.s - s_support
        s_a ~ ratio * s_b
        0 ~ ratio * flange_a.f + flange_b.f
    ]
    extend(ODESystem(eqs, t, sts, [ratio]; name=name), partial_element)
end
