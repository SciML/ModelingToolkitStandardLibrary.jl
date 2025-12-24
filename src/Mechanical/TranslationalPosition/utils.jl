@connector function Flange(; name, s = nothing, f = nothing)
    vars = @variables begin
        s(t) = s
        f(t) = f, [connect = Flow]
    end
    System(Equation[], t, vars, []; name)
end
Base.@doc """
    Flange(;name)

1-dim. translational flange.

# States:
- `s`: [m] Absolute position of flange
- `f`: [N] Cut force into the flange
""" Flange

@connector function Support(; name, s = nothing, f = nothing)
    vars = @variables begin
        s(t) = s
        f(t) = f, [connect = Flow]
    end
    System(Equation[], t, vars, []; name)
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
@component function PartialCompliant(; name, v_a = nothing, v_b = nothing, s_rel = nothing, f = nothing)
    pars = @parameters begin
    end

    systems = @named begin
        flange_a = Flange()
        flange_b = Flange()
    end

    vars = @variables begin
        v_a(t) = v_a
        v_b(t) = v_b
        s_rel(t) = s_rel
        f(t) = f
    end

    equations = Equation[
        D(flange_a.s) ~ v_a,
        D(flange_b.s) ~ v_b,
        D(s_rel) ~ v_b - v_a,
        flange_b.f ~ +f,
        flange_a.f ~ -f
    ]

    return System(equations, t, vars, pars; name, systems)
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
@component function PartialCompliantWithRelativeStates(; name, delta_s = nothing, f = nothing)
    pars = @parameters begin
    end

    systems = @named begin
        flange_a = Flange()
        flange_b = Flange()
    end

    vars = @variables begin
        delta_s(t) = delta_s
        f(t) = f
    end

    equations = Equation[
        delta_s ~ flange_a.s - flange_b.s,
        flange_a.f ~ +f,
        flange_b.f ~ -f
    ]

    return System(equations, t, vars, pars; name, systems)
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
    @variables s_support(t)
    if use_support
        @named support = Support()
        eqs = [support.s ~ s_support
               support.f ~ -flange.f]
        push!(sys, support)
    else
        eqs = [s_support ~ 0]
    end
    return compose(System(eqs, t, [s_support], []; name = name), sys)
end

"""
    PartialElementaryTwoFlangesAndSupport2(; name, use_support = false)

Partial model for a component with two translational 1-dim. flanges and a support used for textual modeling, i.e., for elementary models

# Parameters:

  - `use_support`: If support flange enabled, otherwise implicitly grounded. By default it is `false`

# States:

  - `s_support`: [m] Absolute position of support flange"
"""
@component function PartialElementaryTwoFlangesAndSupport2(; name, use_support = false)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sys = [flange_a, flange_b]
    @variables s_support(t)
    if use_support
        @named support = Support()
        eqs = [support.s ~ s_support
               support.f ~ -flange_a.f - flange_b.f]
        push!(sys, support)
    else
        eqs = [s_support ~ 0]
    end
    return compose(System(eqs, t, [s_support], []; name = name), sys)
end
