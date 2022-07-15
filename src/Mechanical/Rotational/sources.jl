"""
    Torque(;name) 

Input signal acting as external torque on a flange
"""
function Torque(; name, use_support = false)
    @named partial_element = PartialElementaryOneFlangeAndSupport2(use_support = use_support)
    @unpack flange = partial_element
    @named tau = RealInput() # Accelerating torque acting at flange (= -flange.tau)
    eqs = [flange.tau ~ -tau.u]
    return extend(ODESystem(eqs, t, [], []; name = name, systems = [tau]), partial_element)
end
