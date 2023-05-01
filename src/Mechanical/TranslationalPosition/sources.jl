"""
    Force(;name)

Input signal acting as external force on a flange
"""
@component function Force(; name, use_support = false)
    @named partial_element = PartialElementaryOneFlangeAndSupport2(use_support = use_support)
    @unpack flange = partial_element
    @named input = RealInput() # Accelerating force acting at flange (= -flange.tau)
    eqs = [flange.f ~ -input.u]
    return extend(ODESystem(eqs, t, [], []; name = name, systems = [input]),
                  partial_element)
end
