"""
    Fixed(;name, phi0 = 0.0)

Flange fixed in housing at a given angle.

# Connectors:

  - `flange` [Flange](@ref)

# Parameters:

  - `phi0`: [`rad`] Fixed offset angle of housing
"""
@component function Fixed(; phi0 = 0.0, name)
    pars = @parameters begin
        phi0 = phi0, [description = "Fixed offset angle of flange"]
    end

    systems = @named begin
        flange = Flange()
    end

    vars = @variables begin
    end

    equations = Equation[
        flange.phi ~ phi0,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Inertia(;name, J, phi = nothing, w = nothing, a = nothing)

1D-rotational component with inertia.

# States:

  - `phi`: [`rad`] Absolute rotation angle of component
  - `w`: [`rad/s`] Absolute angular velocity of component (= D(phi))
  - `a`: [`rad/s²`] Absolute angular acceleration of component (= D(w))

# Connectors:

  - `flange_a` [Flange](@ref) Left flange
  - `flange_b` [Flange](@ref) Right flange

# Parameters:

  - `J`: [`kg·m²`] Moment of inertia
"""
@component function Inertia(; J = nothing, phi = nothing, w = nothing, a = nothing, name)
    @symcheck J > 0 || throw(ArgumentError("Expected `J` to be positive"))

    pars = @parameters begin
        J = J, [description = "Moment of inertia"]
    end

    systems = @named begin
        flange_a = Flange()
        flange_b = Flange()
    end

    vars = @variables begin
        phi(t) = phi, [description = "Absolute rotation angle", guess = 0.0]
        w(t) = w, [description = "Absolute angular velocity", guess = 0.0]
        a(t) = a, [description = "Absolute angular acceleration", guess = 0.0]
    end

    equations = Equation[
        phi ~ flange_a.phi,
        phi ~ flange_b.phi,
        D(phi) ~ w,
        D(w) ~ a,
        J * a ~ flange_a.tau + flange_b.tau,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Spring(; name, c, phi_rel0 = 0.0)

Linear 1D rotational spring

# States:

  - `phi_rel(t)`: [`rad`] Relative rotation angle (`flange_b.phi - flange_a.phi`)
  - `tau(t)`: [`N.m`] Torque between flanges (`flange_b.tau`)

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)

# Parameters:

  - `c`: [`N.m/rad`] Spring constant
  - `phi_rel0`: [`rad`] Unstretched spring angle. Defaults to 0.0.
"""
@component function Spring(; c = nothing, phi_rel0 = 0.0, name)
    @symcheck c > 0 || throw(ArgumentError("Expected `c` to be positive"))

    @named partial_comp = PartialCompliant()
    @unpack phi_rel, tau = partial_comp

    pars = @parameters begin
        c = c, [description = "Spring constant"]
        phi_rel0 = phi_rel0, [description = "Unstretched spring angle"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        tau ~ c * (phi_rel - phi_rel0),
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, partial_comp)
end

"""
    Damper(; name, d)

Linear 1D rotational damper

# States:

  - `phi_rel(t)`: [`rad`] Relative rotation angle (= flange_b.phi - flange_a.phi)
  - `w_rel(t)`: [`rad/s`] Relative angular velocity (= D(phi_rel))
  - `a_rel(t)`: [`rad/s²`] Relative angular acceleration (= D(w_rel))
  - `tau(t)`: [`N.m`] Torque between flanges (= flange_b.tau)

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)

# Parameters:

  - `d`: [`N.m.s/rad`] Damping constant
"""
@component function Damper(; d = nothing, name)
    @symcheck d > 0 || throw(ArgumentError("Expected `d` to be positive"))

    @named partial_comp = PartialCompliantWithRelativeStates()
    @unpack w_rel, tau = partial_comp

    pars = @parameters begin
        d = d, [description = "Damping constant"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        tau ~ d * w_rel,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, partial_comp)
end
"""
    SpringDamper(; name, d)

Linear 1D rotational spring and damper

# States:

  - `phi_rel(t)`: [`rad`] Relative rotation angle (= flange_b.phi - flange_a.phi)
  - `w_rel(t)`: [`rad/s`] Relative angular velocity (= D(phi_rel))
  - `a_rel(t)`: [`rad/s²`] Relative angular acceleration (= D(w_rel))
  - `tau(t)`: [`N.m`] Torque between flanges (= flange_b.tau)

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)

# Parameters:

  - `d`: [`N.m.s/rad`] Damping constant
  - `c`: [`N.m/rad`] Spring constant
  - `phi_rel0`: [`rad`] Unstretched spring angle. Defaults to 0.0
"""
@component function SpringDamper(; d = nothing, c = nothing, phi_rel0 = 0.0, tau_c = nothing, tau_d = nothing, name)
    @named partial_comp = PartialCompliantWithRelativeStates()
    @unpack phi_rel, w_rel, tau = partial_comp

    pars = @parameters begin
        d = d, [description = "Damping constant"]
        c = c, [description = "Spring constant"]
        phi_rel0 = phi_rel0, [description = "Unstretched spring angle"]
    end

    systems = @named begin
    end

    vars = @variables begin
        tau_c(t) = tau_c, [description = "Spring torque"]
        tau_d(t) = tau_d, [description = "Damper torque"]
    end

    equations = Equation[
        tau_c ~ c * (phi_rel - phi_rel0),
        tau_d ~ d * w_rel,
        tau ~ tau_c + tau_d,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, partial_comp)
end

"""
    IdealGear(; name, ratio, use_support = false)

Ideal gear without inertia.

This element characterizes any type of gear box which is fixed in the ground and which has one driving shaft and one driven shaft.

# States:

  - `phi_a(t)`: [`rad`] Relative angle between shaft a and the support
  - `phi_b(t)`: [`rad`] Relative angle between shaft b and the support

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)
  - `support` [Support](@ref) if `use_support == true`

# Parameters:

  - `ratio`: Transmission ratio (flange_a.phi/flange_b.phi)
  - `use_support`: If support flange enabled, otherwise implicitly grounded. By default it is `false`
"""
@component function IdealGear(; ratio = nothing, use_support = false, phi_a = nothing, phi_b = nothing, name)
    @named partial_element = PartialElementaryTwoFlangesAndSupport2(; use_support)
    @unpack phi_support, flange_a, flange_b = partial_element

    pars = @parameters begin
        ratio = ratio, [description = "Transmission ratio"]
    end

    systems = @named begin
    end

    vars = @variables begin
        phi_a(t) = phi_a, [description = "Relative angle between shaft a and the support", guess = 0.0]
        phi_b(t) = phi_b, [description = "Relative angle between shaft b and the support", guess = 0.0]
    end

    equations = Equation[
        phi_a ~ flange_a.phi - phi_support,
        phi_b ~ flange_b.phi - phi_support,
        phi_a ~ ratio * phi_b,
        0 ~ ratio * flange_a.tau + flange_b.tau,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, partial_element)
end

"""
    RotationalFriction(; name, f, tau_c, w_brk, tau_brk)

Models rotational friction with Stribeck effect, Coulomb friction and viscous friction between the two flanges.
The friction torque is a function of the relative angular velocity between `flange_a` and `flange_b`.

Friction model: "Armstrong, B. and C.C. de Wit, Friction Modeling and Compensation, The Control Handbook, CRC Press, 1995."

# States:

  - `phi_rel(t)`: [`rad`] Relative rotation angle `(= flange_b.phi - flange_a.phi)`
  - `w_rel(t)`: [`rad/s`] Relative angular velocity `(= D(phi_rel))`
  - `a_rel(t)`: [`rad/s²`] Relative angular acceleration `(= D(w_rel))`
  - `tau(t)`: [`N.m`] Torque between flanges `(= flange_b.tau)`

# Connectors:

  - `flange_a` [Flange](@ref)
  - `flange_b` [Flange](@ref)

# Parameters:

  - `f`: [`N⋅m/(rad/s)`] Viscous friction coefficient
  - `tau_c`: [`N⋅m`] Coulomb friction torque
  - `w_brk`: [`rad/s`] Breakaway friction velocity
  - `tau_brk`: [`N⋅m`] Breakaway friction torque
"""
@component function RotationalFriction(; f = nothing, tau_c = nothing, w_brk = nothing, tau_brk = nothing, name)
    @named partial_comp = PartialCompliantWithRelativeStates()
    @unpack w_rel, tau = partial_comp

    pars = @parameters begin
        f = f, [description = "Viscous friction coefficient"]
        tau_c = tau_c, [description = "Coulomb friction torque"]
        w_brk = w_brk, [description = "Breakaway friction velocity"]
        tau_brk = tau_brk, [description = "Breakaway friction torque"]
    end

    str_scale = sqrt(2 * exp(1)) * (tau_brk - tau_c)
    w_st = w_brk * sqrt(2)
    w_coul = w_brk / 10

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        tau ~
            str_scale * (exp(-(w_rel / w_st)^2) * w_rel / w_st) +
            tau_c * tanh(w_rel / w_coul) + f * w_rel,  # Stribeck friction + Coulomb friction + Viscous friction
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, partial_comp)
end
