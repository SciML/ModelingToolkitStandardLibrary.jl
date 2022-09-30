@connector function Flange(; name)
    vars = @variables begin
        s(t)
        f(t), [connect = Flow]
    end
    ODESystem(Equation[], t, vars, [], name = name, defaults = Dict(f => 0.0))
end
Base.@doc """
    Flange(;name)

1-dim. translational flange.

# States:
- `s`: [m] Absolute position of flange
- `f`: [N] Cut force into the flange
""" Flange

@connector function Support(; name)
    sts = @variables begin
        s(t)
        f(t), [connect = Flow]
    end
    ODESystem(Equation[], t, sts, [], name = name, defaults = Dict(s => 0.0, f => 0.0))
end
Base.@doc """
    Support(;name)

Support/housing 1-dim. translational flange.

# States:
- `s`: [m] Absolute position of the support/housing
- `f`: [N] Cut force into the flange
""" Support

"""
    PartialCompliant(;name, s_rel_start=0.0, f_start=0.0)

Partial model for the compliant connection of two translational 1-dim. flanges.

# Parameters:
- `s_rel_start`: [m] Initial relative distance between the flanges 
- `f_start`: [N] Initial force between flanges

# States:
- `s_rel`: [m] Relative distance (= flange_b.s - flange_a.s)
- `f`: [N] Force between flanges (= flange_b.f)
"""
function PartialCompliant(; name, s_rel_start = 0.0, f_start = 0.0)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sts = @variables begin
        v_a(t) = 0
        v_b(t) = 0
        s_rel(t) = s_rel_start
        f(t) = f_start
    end
    eqs = [D(flange_a.s) ~ v_a
           D(flange_b.s) ~ v_b
           D(s_rel) ~ v_b - v_a
           flange_b.f ~ +f
           flange_a.f ~ -f]
    return compose(ODESystem(eqs, t, sts, []; name = name), flange_a, flange_b)
end

"""
    PartialCompliantWithRelativeStates(;name, s_rel_start=0.0, v_rel_start=0.0, a_rel_start=0.0, f_start=0.0)

Partial model for the compliant connection of two translational 1-dim. flanges.

# Parameters:
- `s_rel_start`: [m] Initial relative distance
- `v_rel_start`: [m/s] Initial relative linear velocity (= der(s_rel))
- `a_rel_start`: [m/s²] Initial relative linear acceleration (= der(v_rel))
- `f_start`: [N] Initial force between flanges

# States:
- `s_rel`: [m] Relative distance (= flange_b.phi - flange_a.phi)
- `v_rel`: [m/s] Relative linear velocity (= der(s_rel))
- `a_rel`: [m/s²] Relative linear acceleration (= der(v_rel))
- `f`: [N] Force between flanges (= flange_b.f)
"""
function PartialCompliantWithRelativeStates(; name, delta_s0 = 0.0)
    @named port_a = Flange()
    @named port_b = Flange()
    sts = @variables begin
        delta_s(t) = delta_s0
        f(t) = 0
    end
    eqs = [delta_s ~ port_a.s - port_b.s
           port_a.f ~ +f
           port_b.f ~ -f]
    return compose(ODESystem(eqs, t, sts, []; name = name), port_a, port_b)
end

"""
    PartialElementaryOneFlangeAndSupport2(;name, use_support=false)

Partial model for a component with one translational 1-dim. shaft flange and a support used for textual modeling, i.e., for elementary models

# Parameters:
- `use_support`: If support flange enabled, otherwise implicitly grounded

# States:
- `s_support`: [m] Absolute position of support flange"
"""
function PartialElementaryOneFlangeAndSupport2(; name, use_support = false)
    @named flange = Flange()
    sys = [flange]
    @variables s_support(t)
    if use_support
        @named support = Support()
        eqs = [support.s ~ s_support
               support.f ~ -flange.f]
        push!(sys, support)
    else
        eqs = [s_support ~ 0]
    end
    return compose(ODESystem(eqs, t, [s_support], []; name = name), sys)
end

"""
    PartialElementaryTwoFlangesAndSupport2(;name, use_support=false)

Partial model for a component with two translational 1-dim. flanges and a support used for textual modeling, i.e., for elementary models

# Parameters:
- `use_support`: If support flange enabled, otherwise implicitly grounded

# States:
- `s_support`: [m] Absolute position of support flange"
"""
function PartialElementaryTwoFlangesAndSupport2(; name, use_support = false)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sys = [flange_a, flange_b]
    @variables s_support(t) = 0.0
    if use_support
        @named support = Support()
        eqs = [support.s ~ s_support
               support.f ~ -flange_a.f - flange_b.f]
        push!(sys, support)
    else
        eqs = [s_support ~ 0]
    end
    return compose(ODESystem(eqs, t, [s_support], []; name = name), sys)
end
