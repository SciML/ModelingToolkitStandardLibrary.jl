"""
Constant magnetomotive force.

Parameters:

  - `V_m`: [A] Magnetic potential difference
"""
function ConstantMagneticPotentialDifference(; name, V_m = 1.0)
    port_p = PositiveMagneticPort()
    port_n = NegativeMagneticPort()
    @parameters V_m = V_m
    @variables Phi(t)
    eqs = [V_m ~ port_p.V_m - port_n.V_m
           Phi ~ port_p.Phi
           0 ~ port_p.Phi + port_n.Phi]
    ODESystem(eqs, t, [Phi], [V_m], systems = [port_p, port_n], name = name)
end

"""
Source of constant magnetic flux.

Parameters:

  - `Phi`: [Wb] Magnetic flux
"""
function ConstantMagneticFlux(; name, Phi = 1.0)
    port_p = PositiveMagneticPort()
    port_n = NegativeMagneticPort()
    @parameters Phi = Phi
    @variables V_m(t)
    eqs = [V_m ~ port_p.V_m - port_n.V_m
           Phi ~ port_p.Phi
           0 ~ port_p.Phi + port_n.Phi]
    ODESystem(eqs, t, [Phi], [V_m], systems = [port_p, port_n], name = name)
end
