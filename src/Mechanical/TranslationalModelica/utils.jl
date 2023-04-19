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
    @named flange = Flange()
    extend(ODESystem(Equation[], t, name = name), flange)
end
Base.@doc """
    Support(;name)

Support/housing 1-dim. translational flange.

# States:
- `s`: [m] Absolute position of the support/housing
- `f`: [N] Cut force into the flange
""" Support

function PartialTwoFlanges(; name)
    @named flange_a = Flange() # (left) driving flange (flange axis directed into cut plane, e. g. from left to right)
    @named flange_b = Flange() # (right) driven flange (flange axis directed out of cut plane)
    compose(ODESystem([], t; name), flange_a, flange_b)
end

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
    @named pt = PartialTwoFlanges()
    @unpack flange_a, flange_b = pt
    @variables s_rel(t)=s_rel_start [
        description = "Relative distance between flanges flange_b.s - flange_a.s",
    ]
    @variables f(t)=f_start [
        description = "Force between flanges (positive in direction of flange axis R)",
    ]

    eqs = [s_rel ~ flange_b.s - flange_a.s
           flange_b.f ~ +f
           flange_a.f ~ -f]
    return extend(ODESystem(eqs, t; name = name), pt)
end

"""
    PartialCompliantWithRelativeStates(;name, s_rel_start=0.0, v_rel_start=0.0, f_start=0.0)

Partial model for the compliant connection of two translational 1-dim. flanges.

# Parameters:

  - `s_rel_start`: [m] Initial relative distance
  - `v_rel_start`: [m/s] Initial relative linear velocity (= der(s_rel))

# States:

  - `s_rel`: [m] Relative distance (= flange_b.phi - flange_a.phi)
  - `v_rel`: [m/s] Relative linear velocity (= der(s_rel))
  - `f`: [N] Force between flanges (= flange_b.f)
"""
function PartialCompliantWithRelativeStates(; name, s_rel_start = 0, v_rel_start = 0,
                                            f_start = 0)
    @named pt = PartialTwoFlanges()
    @unpack flange_a, flange_b = pt
    @variables s_rel(t)=s_rel_start [
        description = "Relative distance between flanges flange_b.s - flange_a.s",
    ]
    @variables v_rel(t)=v_rel_start [description = "Relative linear velocity (= D(s_rel))"]
    @variables f(t)=f_start [description = "Forces between flanges (= flange_b.f)"]

    eqs = [s_rel ~ flange_b.s - flange_a.s
           v_rel ~ D(s_rel)
           flange_b.f ~ f
           flange_a.f ~ -f]
    return extend(ODESystem(eqs, t; name = name), pt)
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
    @variables s_support(t) [description = "Absolute position of support flange"]
    @variables s(t) [
        description = "Distance between flange and support (= flange.s - support.s)",
    ]
    eqs = [s ~ flange.s - s_support]
    if use_support
        @named support = Support()
        push!(eqs, support.f ~ -flange.f)
        compose(ODESystem(eqs, t; name = name), flange, support)
    else
        push!(eqs, s_support ~ 0)
        compose(ODESystem(eqs, t; name = name), flange)
    end
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
    @named flange = Flange()

    @variables s_a(t) [description = "Distance between left flange and support"]
    @variables s_b(t) [description = "Distance between right flange and support"]
    @variables s_support(t) [description = "Absolute position of support flange"]

    eqs = [s_a ~ flange_a.s - s_support
           s_b ~ flange_b.s - s_support]
    if use_support
        @named support = Support()
        push!(eqs, support.f ~ -flange_a.f - flange_b.f)
        compose(ODESystem(eqs, t; name = name), flange, support)
    else
        push!(eqs, s_support ~ 0)
        compose(ODESystem(eqs, t; name = name), flange)
    end
end

function PartialRigid(; name, L = 0, s0 = 0)
    @named ptf = PartialTwoFlanges()
    @unpack flange_a, flange_b = ptf
    @variables s(t)=s0 [
        description = "Absolute position of center of component (s = flange_a.s + L/2 = flange_b.s - L/2)",
    ]
    @parameters L=L [
        description = "Length of component, from left flange to right flange (= flange_b.s - flange_a.s)",
    ]
    eqs = [flange_a.s ~ s - L / 2
           flange_b.s ~ s + L / 2]
    return extend(ODESystem(eqs, t; name = name), ptf)
end
