"""
    Force(; name, use_support = false)

Input signal acting as external force on a flange
"""
@component function Force(; name, use_support = false)
    @named partial_element = PartialElementaryOneFlangeAndSupport2(; use_support)
    @unpack flange = partial_element

    pars = @parameters begin
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

"""
    Position(; name, exact = false, f_crit = 50)

Forced movement of a flange according to a reference position

The input signal `s_ref` defines the reference position in [m]. Flange flange is forced to move relative to the support connector according to this reference motion. According to parameter `exact`, this is done in the following way:

- `exact=true`: The reference position is treated exactly. This is only possible, if the input signal is defined by an analytical function which can be differentiated at least twice. If this prerequisite is fulfilled, the Modelica translator will differentiate the input signal twice in order to compute the reference acceleration of the flange.
- `exact=false`: The reference position is filtered and the second derivative of the filtered curve is used to compute the reference acceleration of the flange. This second derivative is not computed by numerical differentiation but by an appropriate realization of the filter. For filtering, a second order Bessel filter is used. The critical frequency (also called cut-off frequency) of the filter is defined via parameter `f_crit` in [Hz]. This value should be selected in such a way that it is higher as the essential low frequencies in the signal.

The input signal can be provided from one of the signal generator blocks of the block library `Blocks.Sources`.
"""
@component function Position(; name, exact = false, f_crit = 50, v = nothing, a = nothing)
    @named ptf = PartialElementaryOneFlangeAndSupport2()
    @unpack s = ptf

    pars = @parameters begin
        f_crit = f_crit
    end

    w_crit = 2π * f_crit
    af = 1.3617
    bf = 0.618

    systems = @named begin
        s_ref = RealInput()
    end

    vars = @variables begin
        v(t) = v
        a(t) = a
    end

    equations = if exact
        Equation[
            s ~ s_ref.u,
            v ~ D(s),
            a ~ D(v),
        ]
    else
        Equation[
            a ~ ((s_ref.u - s) * w_crit - af * v) * (w_crit / bf),
            v ~ D(s),
            a ~ D(v),
        ]
    end

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, ptf)
end
