@connector function MagneticPort(; name, V_m_start = 0.0, Phi_start = 0.0)
    @variables V_m(t) = V_m_start # [Wb] Magnetic potential at the port
    @variables Phi(t)=Phi_start [connect = Flow] # [A] Magnetic flux flowing into the port"
    ODESystem(Equation[], t, [V_m, Phi], []; name = name)
end
Base.@doc "Port for a Magnetic system." MagneticPort

"""
Positive magnetic port
"""
const PositiveMagneticPort = MagneticPort

"""
Negative magnetic port
"""
const NegativeMagneticPort = MagneticPort

"""
    TwoPort(;name, V_m_start=0.0, Phi_start=0.0)

Partial component with magnetic potential difference between two magnetic ports p and n and magnetic flux Phi from p to n.

# Parameters:

  - `V_m_start`: Initial magnetic potential difference between both ports
  - `Phi_start`: Initial magnetic flux from port_p to port_n
"""
@component function TwoPort(; name, V_m_start = 0.0, Phi_start = 0.0)
    @named port_p = PositiveMagneticPort()
    @named port_n = NegativeMagneticPort()
    @variables V_m(t)=V_m_start Phi(t)=Phi_start
    eqs = [V_m ~ port_p.V_m - port_n.V_m
           Phi ~ port_p.Phi
           0 ~ port_p.Phi + port_n.Phi]
    ODESystem(eqs, t, [V_m, Phi], [], systems = [port_p, port_n], name = name)
end
