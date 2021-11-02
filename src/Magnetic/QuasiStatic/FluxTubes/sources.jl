function ConstantMagneticPotentialDifference(;name, V_m=1.0)
    val = V_m
    @named two_port_elementary = TwoPortElementary()
    @unpack port_p, port_n = two_port_elementary
    @parameters V_m
    @variables Phi(t)
    eqs = [
        V_m ~ port_p.V_m - port_n.V_m,
        Phi ~ port_p.Phi,
        0 ~ port_p.Phi + port_n.Phi,
    ]
    extend(ODESystem(eqs, t, [Phi], [V_m], systems=[port_p, port_n], defaults=Dict(V_m => val), name=name), two_port_elementary)
end

function ConstantMagneticFlux(;name, Phi=1.0)
    val = Phi
    @named two_port_elementary = TwoPortElementary()
    @unpack port_p, port_n = two_port_elementary
    @parameters Phi
    @variables V_m(t)
    eqs = [
        V_m ~ port_p.V_m - port_n.V_m,
        Phi ~ port_p.Phi,
        0 ~ port_p.Phi + port_n.Phi,
    ]
    extend(ODESystem(eqs, t, [V_m], [Phi], systems=[port_p, port_n], defaults=Dict(Phi => val), name=name), two_port_elementary)
end
