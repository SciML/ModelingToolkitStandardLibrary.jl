"""
    Force(;name)

Input signal acting as external force on a flange
"""
@component function Force(; name, use_support = false)
    @named partial_element = PartialElementaryOneFlangeAndSupport2(use_support = use_support)
    @unpack flange = partial_element
    @named f = RealInput() # Accelerating force acting at flange (= -flange.tau)
    eqs = [flange.f ~ -f.u]
    return extend(ODESystem(eqs, t, [], []; name = name, systems = [f]), partial_element)
end
