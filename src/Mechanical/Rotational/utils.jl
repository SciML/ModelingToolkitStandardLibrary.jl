@connector function Flange(; name)
    sts = @variables(phi(t), [description = "Rotation angle of flange $name"],
                     tau(t), [connect = Flow, description = "Cut torque in flange $name"],)
    ODESystem(Equation[], t, sts, [], name = name, defaults = Dict(phi => 0.0, tau => 0.0))
end
Base.@doc """
    Support(;name)

1-dim. rotational flange of a shaft.

# States:
- `phi(t)`: [`rad`] Absolute rotation angle of flange
- `tau(t)`: [`N.m`] Cut torque in the flange
""" Flange

@connector function Support(; name)
    @named flange = Flange()
    extend(ODESystem(Equation[], t, [], [], name = name), flange)
end


Base.@doc """
    InternalSupport(;name, tau)

1-dim. rotational flange of a shaft.

- `tau`: External support torque (must be computed via torque balance in model where InternalSupport is used; = flange.tau)

# States:
- `phi(t)`: [`rad`] Absolute rotation angle of flange
- `tau(t)`: [`N.m`] Cut torque in the flange
""" Flange

@connector function InternalSupport(; name, tau)
    @named flange = Flange()
    @variables phi(t)=0 [description = "Rotation angle of support $name"]
    # tau(t), [connect = Flow, description = "Cut torque in support $name"],)
    equations = [
        flange.tau ~ tau
        flange.phi ~ phi
    ]
    ODESystem(equations, t, name = name, systems=[flange])
end

Base.@doc """
    Support(;name)

Support/housing of a 1-dim. rotational shaft

# States:
- `phi(t)`: [`rad`] Absolute rotation angle of the support/housing
- `tau(t)`: [`N.m`] Cut torque in the support/housing
""" Support

"""
    PartialCompliant(;name, phi_rel_start=0.0, tau_start=0.0)

Partial model for the compliant connection of two rotational 1-dim. shaft flanges.

# States:

  - `phi_rel(t)`: [`rad`] Relative rotation angle (`flange_b.phi - flange_a.phi`)
  - `tau(t)`: [`N.m`] Torque between flanges (`flange_b.tau`)

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)

# Parameters:

  - `phi_rel_start`: [`rad`] Initial relative rotation angle
  - `tau_start`: [`N.m`] Initial torque between flanges
"""
function PartialCompliant(; name, phi_rel_start = 0.0, tau_start = 0.0)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sts = @variables(phi_rel(t)=phi_rel_start,
                     [description = "Relative rotation angle between flanges"],
                     tau(t)=tau_start, [description = "Torque between flanges"])
    eqs = [phi_rel ~ flange_b.phi - flange_a.phi
           flange_b.tau ~ tau
           flange_a.tau ~ -tau]
    return compose(ODESystem(eqs, t, sts, []; name = name), flange_a, flange_b)
end

"""
    PartialCompliantWithRelativeStates(;name, phi_rel_start=0.0, tau_start=0.0)

Partial model for the compliant connection of two rotational 1-dim. shaft flanges where the relative angle and speed are used as preferred states

# States:

  - `phi_rel(t)`: [`rad`] Relative rotation angle (= flange_b.phi - flange_a.phi)
  - `w_rel(t)`: [`rad/s`] Relative angular velocity (= D(phi_rel))
  - `a_rel(t)`: [`rad/s²`] Relative angular acceleration (= D(w_rel))
  - `tau(t)`: [`N.m`] Torque between flanges (= flange_b.tau)

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)

# Parameters:

  - `phi_rel_start`: [`rad`] Initial relative rotation angle
  - `w_rel_start`: [`rad/s`] Initial relative angular velocity (= D(phi_rel))
  - `a_rel_start`: [`rad/s²`] Initial relative angular acceleration (= D(w_rel))
  - `tau_start`: [`N.m`] Initial torque between flanges
"""
function PartialCompliantWithRelativeStates(; name, phi_rel_start = 0.0, w_start = 0.0,
                                            a_start = 0.0, tau_start = 0.0)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sts = @variables(phi_rel(t)=phi_rel_start,
                     [description = "Relative rotation angle between flanges"],
                     w_rel(t)=w_start,
                     [description = "Relative angular velocity between flanges"],
                     a_rel(t)=a_start,
                     [description = "Relative angular acceleration between flanges"],
                     tau(t)=tau_start, [description = "Torque between flanges"],)
    eqs = [phi_rel ~ flange_b.phi - flange_a.phi
           D(phi_rel) ~ w_rel
           D(w_rel) ~ a_rel
           flange_b.tau ~ tau
           flange_a.tau ~ -tau]
    return compose(ODESystem(eqs, t, sts, []; name = name), flange_a, flange_b)
end

"""
    PartialElementaryOneFlangeAndSupport2(;name, use_support=false)

Partial model for a component with one rotational 1-dim. shaft flange and a support used for textual modeling, i.e., for elementary models

# States:

  - `phi_support(t)`: [`rad`] Absolute angle of support flange"

# Connectors:

  - `flange` [Flange](@ref)

# Parameters:

  - `use_support`: If support flange enabled, otherwise implicitly grounded
"""
function PartialElementaryOneFlangeAndSupport2(; name, use_support = false)
    @named flange = Flange()
    sys = [flange]
    @variables phi_support(t)=0.0 [description = "Absolute angle of support flange"]
    if use_support
        @named support = Support()
        eqs = [support.phi ~ phi_support
               support.tau ~ -flange.tau]
        push!(sys, support)
    else
        eqs = [phi_support ~ 0]
    end
    return compose(ODESystem(eqs, t, [phi_support], []; name = name), sys)
end

"""
    PartialElementaryTwoFlangesAndSupport2(;name, use_support=false)

Partial model for a component with two rotational 1-dim. shaft flanges and a support used for textual modeling, i.e., for elementary models

# States:

  - `phi_support(t)`: [`rad`] Absolute angle of support flange

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)
  - `support` [Support](@ref)  if `use_support == true`

# Parameters:

  - `use_support`: If support flange enabled, otherwise implicitly grounded
"""
function PartialElementaryTwoFlangesAndSupport2(; name, use_support = false)
    @named flange_a = Flange()
    @named flange_b = Flange()
    sys = [flange_a, flange_b]
    @variables phi_support(t)=0.0 [description = "Absolute angle of support flange"]
    if use_support
        @named support = Support()
        eqs = [support.phi ~ phi_support
               support.tau ~ -flange_a.tau - flange_b.tau]
        push!(sys, support)
    else
        eqs = [phi_support ~ 0]
    end
    return compose(ODESystem(eqs, t, [phi_support], []; name = name), sys)
end
