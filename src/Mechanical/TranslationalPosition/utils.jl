@connector Flange begin
    s(t), [description = "Absolute position of flange", unit = u"m"]
    f(t), [connect = Flow, description = " Cut force into the flange", unit = u"N"]
end
Base.@doc """
    Flange(;name)

1-dim. translational flange.

# States:
- `s`: [m] Absolute position of flange
- `f`: [N] Cut force into the flange
""" Flange

@connector Support begin
    s(t), [description = "Absolute position of the support/housing", unit = u"m"]
    f(t), [connect = Flow, description = "Cut force into the flange", unit = u"N"]
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

# States:

  - `s_rel`: [m] Relative distance (= flange_b.s - flange_a.s). It accepts an initial value, which defaults to 0.0.
  - `f`: [N] Force between flanges (= flange_b.f). It accepts an initial value, which defaults to 0.0.
"""
@mtkmodel PartialCompliant begin
    @components begin
        flange_a = Flange()
        flange_b = Flange()
    end
    @variables begin
        v_a(t) = 0.0, [description = "Velocity", unit = u"m/s"]
        v_b(t) = 0.0, [description = "Velocity", unit = u"m/s"]
        s_rel(t) = 0.0, [description = "Relative distance ", unit = u"m"]
        f(t) = 0.0, [description = "Force between flanges", unit = u"N"]
    end
    @equations begin
        D(flange_a.s) ~ v_a
        D(flange_b.s) ~ v_b
        D(s_rel) ~ v_b - v_a
        flange_b.f ~ +f
        flange_a.f ~ -f
    end
end

"""
    PartialCompliantWithRelativeStates(;name, s_rel_start=0.0, v_rel_start=0.0, a_rel_start=0.0, f_start=0.0)

Partial model for the compliant connection of two translational 1-dim. flanges.

  # States:

  - `delta_s`: [m]
  - `f`: [N] Force between flanges (= flange_b.f)
"""
@mtkmodel PartialCompliantWithRelativeStates begin
    @components begin
        flange_a = Flange()
        flange_b = Flange()
    end
    @variables begin
        delta_s(t) = 0.0, [unit = u"m"]
        f(t) = 0.0, [description = "Force between flanges", unit = u"N"]
    end
    @equations begin
        delta_s ~ flange_a.s - flange_b.s
        flange_a.f ~ +f
        flange_b.f ~ -f
    end
end

"""
    PartialElementaryOneFlangeAndSupport2(;name, use_support=false)

Partial model for a component with one translational 1-dim. shaft flange and a support used for textual modeling, i.e., for elementary models

# Parameters:

  - `use_support`: If support flange enabled, otherwise implicitly grounded

# States:

  - `s_support`: [m] Absolute position of support flange"
"""
@component function PartialElementaryOneFlangeAndSupport2(; name, use_support = false)
    @named flange = Flange()
    sys = [flange]
    @variables s_support(t),
    [description = "Absolute position of support flange", unit = u"m"]
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
    PartialElementaryTwoFlangesAndSupport2(; name, use_support = false)

Partial model for a component with two translational 1-dim. flanges and a support used for textual modeling, i.e., for elementary models

# Parameters:

  - `use_support`: If support flange enabled, otherwise implicitly grounded

# States:

  - `s_support`: [m] Absolute position of support flange"
"""
@component function PartialElementaryTwoFlangesAndSupport2(; name, use_support = false)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sys = [flange_a, flange_b]
    @variables s_support(t) = 0.0 [description = "Absolute position of support flange", unit = u"m"]
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
