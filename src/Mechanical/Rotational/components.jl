"""
    Fixed(;name, phi0=0.0)

Flange fixed in housing at a given angle.

# Connectors:
- `flange` [Flange](@ref)

# Parameters:
- `phi0`: [`rad`] Fixed offset angle of housing
"""
function Fixed(;name, phi0=0.0)
    @named flange = Flange()
    @parameters phi0=phi0 
    eqs = [flange.phi ~ phi0]
    return compose(ODESystem(eqs, t, [], [phi0]; name=name), flange)
end

"""
    Inertia(;name, J, phi_start=0.0, w_start=0.0, a_start=0.0)

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
- `phi_start`: [`rad`] Initial value of absolute rotation angle of component 
- `w_start`: [`rad/s`] Initial value of absolute angular velocity of component
- `a_start`: [`rad/s²`] Initial value of absolute angular acceleration of component
"""
function Inertia(;name, J, phi_start=0.0, w_start=0.0, a_start=0.0)
    @named flange_a = Flange()
    @named flange_b = Flange()
    @parameters J=J
    sts = @variables begin
        phi(t)=phi_start
        w(t)=w_start
        a(t)=a_start
    end
    eqs = [
        phi ~ flange_a.phi
        phi ~ flange_b.phi
        D(phi) ~ w 
        D(w) ~ a 
        J*a ~ flange_a.tau + flange_b.tau
    ]
    return compose(ODESystem(eqs, t, sts, [J]; name=name), flange_a, flange_b)
end

"""
    Spring(;name, c, phi_rel0=0.0)

Linear 1D rotational spring

# States:
- `phi_rel(t)`: [`rad`] Relative rotation angle (`flange_b.phi - flange_a.phi`)
- `tau(t)`: [`N.m`] Torque between flanges (`flange_b.tau`)

# Connectors:
- `flange_a` [Flange](@ref)
- `flange_b` [Flange](@ref)

# Parameters:
- `c`: [`N.m/rad`] Spring constant
- `phi_rel0`: [`rad`] Unstretched spring angle
"""
function Spring(;name, c, phi_rel0=0.0)
    @named partial_comp = PartialCompliant()
    @unpack phi_rel, tau = partial_comp
    pars = @parameters begin
        c=c 
        phi_rel0=phi_rel0  
    end
    eqs = [tau ~ c*(phi_rel - phi_rel0)]
    extend(ODESystem(eqs, t, [], pars; name=name), partial_comp)
end

"""
    Damper(;name, d) 

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
function Damper(;name, d) 
    @named partial_comp = PartialCompliantWithRelativeStates()
    @unpack w_rel, tau = partial_comp
    pars = @parameters d=d 
    eqs = [tau ~ d*w_rel]
    extend(ODESystem(eqs, t, [], pars; name=name), partial_comp)
end

"""
    IdealGear(;name, ratio, use_support=false) 

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
- `use_support`: If support flange enabled, otherwise implicitly grounded
"""
function IdealGear(;name, ratio, use_support=false) 
    @named partial_element = PartialElementaryTwoFlangesAndSupport2(use_support=use_support)
    @unpack phi_support, flange_a, flange_b = partial_element
    @parameters ratio=ratio
    sts = @variables phi_a(t)=0.0 phi_b(t)=0.0
    eqs = [ 
        phi_a ~ flange_a.phi - phi_support
        phi_b ~ flange_b.phi - phi_support
        phi_a ~ ratio*phi_b
        0 ~ ratio*flange_a.tau + flange_b.tau
    ]
    extend(ODESystem(eqs, t, sts, [ratio]; name=name), partial_element)
end


"""
    RotationalFriction(;name, f, tau_c, w_brk, tau_brk) 

Models rotational friction with Stribeck effect, Coulomb friction and viscous friction between the two flanges.
The friction torque is a function of the relative angular velocity between flange_a and flange_b.

Friction model: "Armstrong, B. and C.C. de Wit, Friction Modeling and Compensation, The Control Handbook, CRC Press, 1995."

# States:
- `phi_rel(t)`: [`rad`] Relative rotation angle (= flange_b.phi - flange_a.phi)
- `w_rel(t)`: [`rad/s`] Relative angular velocity (= D(phi_rel))
- `a_rel(t)`: [`rad/s²`] Relative angular acceleration (= D(w_rel))
- `tau(t)`: [`N.m`] Torque between flanges (= flange_b.tau)

# Connectors:
- `flange_a` [Flange](@ref)
- `flange_b` [Flange](@ref)

# Parameters:
- `f`: [`N⋅m/(rad/s)`] Viscous friction coefficient 
- `tau_c`: [`N⋅m`] Coulomb friction torque
- `w_brk`: [`rad/s`] Breakaway friction velocity 
- `tau_brk`: [`N⋅m`] Breakaway friction torque
"""
function RotationalFriction(;name, f, tau_c, w_brk, tau_brk) 
    @named partial_comp = PartialCompliantWithRelativeStates()
    @unpack w_rel, tau = partial_comp
    pars = @parameters f=f tau_c=tau_c w_brk=w_brk tau_brk=tau_brk

    str_scale = sqrt(2*exp(1)) * (tau_brk - tau_c)
    w_st = w_brk * sqrt(2)
    w_coul = w_brk / 10                     

    eqs = [ 
        tau ~ str_scale * (exp(-(w_rel/w_st)^2) * w_rel / w_st) + tau_c * tanh(w_rel / w_coul) + f * w_rel # Stribeck friction + Coulomb friction + Viscous friction
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), partial_comp)
end
