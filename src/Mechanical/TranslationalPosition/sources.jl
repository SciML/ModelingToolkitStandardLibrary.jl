"""
    Force(; name, use_support = false)

Input signal acting as external force on a flange
"""
@component function Force(; name, use_support = false, s = 0)
    @named partial_element = PartialElementaryOneFlangeAndSupport2(; use_support)
    @unpack flange = partial_element

    pars = @parameters begin
        s = s
    end

    systems = @named begin
        f = RealInput()
    end

    vars = @variables begin
    end

    equations = Equation[
        flange.f ~ -f.u,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, partial_element)
end
