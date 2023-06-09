
function PartialTorque(; name, use_support = false)
    @named partial_element = PartialElementaryOneFlangeAndSupport2(use_support = use_support)
    @unpack flange, phi_support = partial_element
    @variables phi(t) [
        description = "Angle of flange with respect to support (= flange.phi - support.phi)",
    ]
    eqs = [phi ~ flange.phi - phi_support]
    return extend(ODESystem(eqs, t; name = name), partial_element)
end

"""
    Torque(; name, use_support=false)

Input signal acting as external torque on a flange

# States:

  - `phi_support(t)`: [`rad`] Absolute angle of support flange

# Connectors:

  - `flange` [Flange](@ref)
  - `tau` [RealInput](@ref)  Accelerating torque acting at flange `-flange.tau`

# Parameters:

  - `use_support`
"""
@component function Torque(; name, use_support = false)
    @named partial_element = PartialElementaryOneFlangeAndSupport2(use_support = use_support)
    @unpack flange = partial_element
    @named tau = RealInput()
    eqs = [flange.tau ~ -tau.u]
    return extend(ODESystem(eqs, t, [], []; name = name, systems = [tau]), partial_element)
end

"""
    ConstantTorque(; name, tau_constant, use_support = false)

Constant torque source

# State variables:
  
- `phi_support(t)`: [`rad`] Absolute angle of support flange, only available if `use_support = true`
- `tau`: Accelerating torque acting at flange (= -flange.tau)
- `w`: Angular velocity of flange with respect to support (= der(phi))

# Connectors:
- `flange` [Flange](@ref)

# Arguments:
- `tau_constant`: The constant torque applied by the source
- `use_support`: Whether or not an internal support flange is added, defaults to false.
"""
function ConstantTorque(; name, tau_constant, use_support = false)
    @named partial_element = PartialTorque(; use_support)
    @unpack flange, phi = partial_element
    @parameters tau_constant=tau_constant [
        description = "Constant torque (if negative, torque is acting as load in positive direction of rotation)",
    ]
    @variables tau(t) [description = "Accelerating torque acting at flange (= -flange.tau)"]
    @variables w(t) [
        description = "Angular velocity of flange with respect to support (= der(phi))",
    ]
    eqs = [w ~ D(phi)
        tau ~ -flange.tau
        tau ~ tau_constant]
    return extend(ODESystem(eqs, t; name = name), partial_element)
end

"""
    Speed(; name, use_support=false, exact=false, f_crit=50) 

Forced movement of a flange according to a reference angular velocity signal

# States:

  - `phi_support(t)`: [`rad`] Absolute angle of support flange"

# Connectors:

  - `flange` [Flange](@ref)
  - `w_ref` [RealInput](@ref) Reference angular velocity of flange with respect to support as input signal needs to be continuously differential

# Parameters:

  - `use_support`: If support flange enabled, otherwise implicitly grounded
  - `exact`: true/false exact treatment/filtering the input signal
  - `tau_filt`: [`rad/s`] if exact=false, Time constant of low-pass filter to filter input signal
"""
@component function Speed(; name, use_support = false, exact = false, tau_filt = 50)
    @named partial_element = PartialElementaryOneFlangeAndSupport2(use_support = use_support)
    @unpack flange, phi_support = partial_element
    @named w_ref = RealInput()
    @variables phi(t)=0.0 w(t)=0.0 a(t)=0.0
    eqs = [phi ~ flange.phi - phi_support
        D(phi) ~ w]
    if exact
        pars = []
        push!(eqs, w ~ w_ref.u)
        push!(eqs, a ~ 0)
    else
        pars = @parameters tau_filt = tau_filt
        push!(eqs, D(w) ~ a)
        push!(eqs, a ~ (w_ref.u - w) * tau_filt)
    end
    return extend(ODESystem(eqs, t, [phi, w, a], pars; name = name, systems = [w_ref]),
        partial_element)
end
