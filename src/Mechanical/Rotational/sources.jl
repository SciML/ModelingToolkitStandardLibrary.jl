"""
    Torque(;name) 

Input signal acting as external torque on a flange

# States:
- `phi_support(t)`: [`rad`] Absolute angle of support flange"

# Connectors:
- `flange` [Flange](@ref) 
- `tau` [RealInput](@ref)  Accelerating torque acting at flange `-flange.tau`

# Parameters:
- `use_support`
"""
function Torque(;name, use_support=false) 
    @named partial_element = PartialElementaryOneFlangeAndSupport2(use_support=use_support)
    @unpack flange = partial_element
    @named tau = RealInput()
    eqs = [flange.tau ~ -tau.u]
    return extend(ODESystem(eqs, t, [], []; name=name, systems=[tau]), partial_element)
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
function Speed(;name, use_support=false, exact=false, tau_filt=50) 
    @named partial_element = PartialElementaryOneFlangeAndSupport2(use_support=use_support)
    @unpack flange, phi_support = partial_element
    @named w_ref = RealInput()
    @variables phi(t)=0.0 w(t)=0.0 a(t)=0.0
    eqs = [
        phi ~ flange.phi - phi_support
        D(phi) ~ w
    ]
    if exact
        pars = []
        push!(eqs, w ~ w_ref.u)
        push!(eqs, a ~ 0)
    else
        pars = @parameters tau_filt=tau_filt
        push!(eqs, D(w) ~ a)
        push!(eqs, a ~ (w_ref.u - w) * tau_filt)
    end
    return extend(ODESystem(eqs, t, [phi, w, a], pars; name=name, systems=[w_ref]), partial_element)
end