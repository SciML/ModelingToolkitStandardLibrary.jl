"""
    Fixed(;name, phi0=0.0)

Flange fixed in housing at a given angle.

# Connectors:
- `flange` flange

# Parameters:
- `phi0`: [rad] Fixed offset angle of housing
"""
function Fixed(;name, phi0=0.0)
    @named flange = Flange()
    @parameters phi0=phi0 
    eqs = [flange.phi ~ phi0]
    return compose(ODESystem(eqs, t, [], [phi0]; name=name), flange)
end

"""
    Inertia(;name, J=1.0, phi_start=0.0, w_start=0.0, a_start=0.0)

1D-rotational component with inertia.

# States: 
- `phi`: [rad] Absolute rotation angle of component 
- `w`: [rad/s] Absolute angular velocity of component (= der(phi)) 
- `a`: [rad/s²] Absolute angular acceleration of component (= der(w)) 

# Connectors:
- `flange_a` Left flange
- `flange_b` Right flange

# Parameters:
- `J`: [kg·m²] Moment of inertia 
- `phi_start`: [rad] Initial value of absolute rotation angle of component 
- `w_start`: [rad/s] Initial value of absolute angular velocity of component
- `a_start`: [rad/s²] Initial value of absolute angular acceleration of component
"""
function Inertia(;name, J=1.0, phi_start=0.0, w_start=0.0, a_start=0.0)
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
- See [PartialCompliant](@ref)

# Connectors:
- See [PartialCompliant](@ref)

# Parameters:
- `c`: [N.m/rad] Spring constant
- `phi_rel0`: Unstretched spring angle
"""
function Spring(;name, c=1.0e5, phi_rel0=0.0)
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
    Damper(;name, d=0.0) 

Linear 1D rotational damper

# States:
- See [PartialCompliantWithRelativeStates](@ref)

# Connectors:
- See [PartialCompliantWithRelativeStates](@ref)

# Parameters:
- `d`: [N.m.s/rad] Damping constant
"""
function Damper(;name, d=0.0) 
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
    RotationalFriction(;name, f=0.001, tau_c=20, w_brk=0.1, tau_brk=25, use_support=false) 

Models rotational friction with Stribeck friction, Coulomb friction and viscous friction.

Friction model: "Armstrong, B. and C.C. de Wit, Friction Modeling and Compensation, The Control Handbook, CRC Press, 1995."

# States:
- `phi(t)`: [`rad`] 
- `w(t)`: [`rad/s`]
- `tau(t)`: [`N⋅m`] Friction torque

# Connectors:
- See [PartialElementaryTwoFlangesAndSupport2](@ref)

# Parameters:
- `f`: [`N⋅m/(rad/s)`] Viscous friction coefficient 
- `tau_c`: [`N⋅m`] Coulomb friction torque
- `w_brk`: [`rad/s`] Breakaway friction velocity 
- `tau_brk`: [`N⋅m`] Breakaway friction torque
"""
function RotationalFriction(;name, f=0.001, tau_c=20, w_brk=0.1, tau_brk=25, use_support=false) 
    @named partial_element = PartialElementaryTwoFlangesAndSupport2(use_support=use_support)
    @unpack phi_support, flange_a, flange_b = partial_element
    pars = @parameters f=f tau_c=tau_c w_brk=w_brk tau_brk=tau_brk
    sts = @variables phi(t)=0.0 w(t)=0.0 tau(t)=0.0

    w_st = w_brk * sqrt(2)
    w_coul = w_brk / 10

    eqs = [ 
        D(phi) ~ w
        phi ~ flange_a.phi - phi_support
        flange_b.phi ~ flange_a.phi
        flange_a.tau + flange_b.tau - tau ~ 0
        tau ~ sqrt(2*exp(1)) * (tau_brk - tau_c) * exp(-(w/w_st)^2) * w / w_st + tau_c * tanh(w/w_coul) + f * w # Stribeck friction + Coulomb friction + Viscous friction
    ]
    extend(ODESystem(eqs, t, sts, pars; name=name), partial_element)
end