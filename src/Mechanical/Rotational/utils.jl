@connector function Flange(;name)
    sts = @variables begin
        phi(t)
        tau(t), [connect=Flow]
    end 
    ODESystem(Equation[], t, sts, [], name=name, defaults=Dict(phi=>0.0, tau=>0.0))
end
Base.@doc """
    Flange(;name)

1-dim. rotational flange of a shaft.

# States:
- `phi`: [rad] Absolute rotation angle of flange
- `tau`: [N.m] Cut torque in the flange
""" Flange

@connector function Support(;name)
    sts = @variables begin
        phi(t)
        tau(t), [connect=Flow]
    end 
    ODESystem(Equation[], t, sts, [], name=name, defaults=Dict(phi=>0.0, tau=>0.0))
end
Base.@doc """
    Support(;name)

Support/housing of a 1-dim. rotational shaft

# States:
- `phi`: [rad] Absolute rotation angle of the support/housing
- `tau`: [N.m] Cut torque in the support/housing
""" Support

"""
    PartialCompliant(;name, phi_rel_start=0.0, tau_start=0.0)

Partial model for the compliant connection of two rotational 1-dim. shaft flanges.

# Parameters:
- `phi_rel_start`: [rad] Initial relative rotation angle
- `tau_start`: [N.m] Initial torque between flanges

# States:
- `phi_rel`: [rad] Relative rotation angle (= flange_b.phi - flange_a.phi)
- `tau`: [N.m] Torque between flanges (= flange_b.tau)
"""
function PartialCompliant(;name, phi_rel_start=0.0, tau_start=0.0)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sts = @variables begin
        phi_rel(t)=phi_rel_start
        tau(t)=tau_start
    end
    eqs = [
        phi_rel ~ flange_b.phi - flange_a.phi
        flange_b.tau ~ tau
        flange_a.tau ~ -tau
    ]
    return compose(ODESystem(eqs, t, sts, []; name=name), flange_a, flange_b)
end

"""
    PartialCompliantWithRelativeStates(;name, phi_rel_start=0.0, tau_start=0.0)

Partial model for the compliant connection of two rotational 1-dim. shaft flanges where the relative angle and speed are used as preferred states

# Parameters:
- `phi_rel_start`: [rad] Initial relative rotation angle
- `w_rel_start`: [rad/s] Initial relative angular velocity (= der(phi_rel))
- `a_rel_start`: [rad/s²] Initial relative angular acceleration (= der(w_rel))
- `tau_start`: [N.m] Initial torque between flanges

# States:
- `phi_rel`: [rad] Relative rotation angle (= flange_b.phi - flange_a.phi)
- `w_rel`: [rad/s] Relative angular velocity (= der(phi_rel))
- `a_rel`: [rad/s²] Relative angular acceleration (= der(w_rel))
- `tau`: [N.m] Torque between flanges (= flange_b.tau)
"""
function PartialCompliantWithRelativeStates(;name, phi_rel_start=0.0, w_start=0.0, a_start=0.0, tau_start=0.0)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sts = @variables begin
        phi_rel(t)=phi_rel_start
        w_rel(t)=w_start
        a_rel(t)=a_start
        tau(t)=tau_start
    end
    eqs = [
        phi_rel ~ flange_b.phi - flange_a.phi
        D(phi_rel) ~ w_rel
        D(w_rel) ~ a_rel
        flange_b.tau ~ tau
        flange_a.tau ~ -tau
    ]
    return compose(ODESystem(eqs, t, sts, []; name=name), flange_a, flange_b)
end

"""
    PartialElementaryOneFlangeAndSupport2(;name, use_support=false)

Partial model for a component with one rotational 1-dim. shaft flange and a support used for textual modeling, i.e., for elementary models

# Parameters:
- `use_support`: If support flange enabled, otherwise implicitly grounded

# States:
- `phi_support`: [rad] Absolute angle of support flange"
"""
function PartialElementaryOneFlangeAndSupport2(;name, use_support=false)
    @named flange = Flange()
    sys = [flange]
    @variables phi_support(t)
    if use_support
        @named support = Support() 
        eqs = [
            support.phi ~ phi_support
            support.tau ~ -flange.tau
        ]
        push!(sys, support)
    else
        eqs = [phi_support ~ 0]
    end
    return compose(ODESystem(eqs, t, [phi_support], []; name=name), sys)
end

"""
    PartialElementaryTwoFlangesAndSupport2(;name, use_support=false)

Partial model for a component with two rotational 1-dim. shaft flanges and a support used for textual modeling, i.e., for elementary models

# Parameters:
- `use_support`: If support flange enabled, otherwise implicitly grounded

# States:
- `phi_support`: [rad] Absolute angle of support flange"
"""
function PartialElementaryTwoFlangesAndSupport2(;name, use_support=false)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sys = [flange_a, flange_b]
    @variables phi_support(t)=0.0
    if use_support
        @named support = Support() 
        eqs = [
            support.phi ~ phi_support
            support.tau ~ -flange_a.tau - flange_b.tau
        ]
        push!(sys, support)
    else
        eqs = [phi_support ~ 0]
    end
    return compose(ODESystem(eqs, t, [phi_support], []; name=name), sys)
end