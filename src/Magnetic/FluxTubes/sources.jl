function ConstantMagneticPotentialDifference(;name, 
    V_m=1.0,
    )   
    @named twoport = TwoPortElementary()
    @unpack v, i = twoport
    pars = @parameters V_m=V_m
    eqs = [
        V_m,
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), twoport)
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
