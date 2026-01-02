"""
    ConstantMagneticPotentialDifference(; name, V_m = 0.0)

Constant magnetomotive force.

Parameters:

  - `V_m`: [A] Magnetic potential difference
"""
@component function ConstantMagneticPotentialDifference(; name, V_m = 0.0, Phi = nothing)
    pars = @parameters begin
        V_m = V_m, [description = "Magnetic potential difference"]
    end

    systems = @named begin
        port_p = PositiveMagneticPort()
        port_n = NegativeMagneticPort()
    end

    vars = @variables begin
        Phi(t) = Phi
    end

    equations = Equation[
        V_m ~ port_p.V_m - port_n.V_m,
        Phi ~ port_p.Phi,
        0 ~ port_p.Phi + port_n.Phi
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    ConstantMagneticFlux(; name, Phi = 0.0)

Source of constant magnetic flux.

Parameters:

  - `Phi`: [Wb] Magnetic flux
"""
@component function ConstantMagneticFlux(; name, Phi = 0.0, V_m = nothing)
    pars = @parameters begin
        Phi = Phi, [description = "Magnetic flux"]
    end

    systems = @named begin
        port_p = PositiveMagneticPort()
        port_n = NegativeMagneticPort()
    end

    vars = @variables begin
        V_m(t) = V_m
    end

    equations = Equation[
        V_m ~ port_p.V_m - port_n.V_m,
        Phi ~ port_p.Phi,
        0 ~ port_p.Phi + port_n.Phi
    ]

    return System(equations, t, vars, pars; name, systems)
end
