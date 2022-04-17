@connector function MagneticPort(;name)
    sts = @variables begin
        V_m(t) # [Wb] Magnetic potential at the port
        Phi(t), [connect=Flow] # [A] Magnetic flux flowing into the port"
    end
    ODESystem(Equation[], t, sts, []; name=name)
end
Base.@doc "Port for a Magnetic system." MagneticPort

const PositiveMagneticPort = MagneticPort
const NegativeMagneticPort = MagneticPort

function TwoPort(;name)
    port_p = PositiveMagneticPort()
    port_n = NegativeMagneticPort()
    @variables V_m(t) Phi(t)
    eqs = [
        V_m ~ port_p.V_m - port_n.V_m
        Phi ~ port_p.Phi
        0 ~ port_p.Phi + port_n.Phi
    ]
    ODESystem(eqs, t, [V_m, Phi], [], systems=[port_p, port_n], name=name)
end