@connector function MagneticPort(; name, V_m = nothing, Phi = nothing)
    vars = @variables begin
        V_m(t) = V_m, [description = "Magnetic potential at the port"]
        Phi(t) = Phi, [connect = Flow, description = "Magnetic flux flowing into the port"]
    end
    System(Equation[], t, vars, []; name)
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
    TwoPort(; name, V_m = 0.0, Phi = 0.0)

Partial component with magnetic potential difference between two magnetic ports p and n and magnetic flux Phi from p to n.

# Parameters:

  - `V_m`: Initial magnetic potential difference between both ports
  - `Phi`: Initial magnetic flux from port_p to port_n
"""
@component function TwoPort(; name, V_m = nothing, Phi = nothing)
    pars = @parameters begin
    end

    systems = @named begin
        port_p = PositiveMagneticPort()
        port_n = NegativeMagneticPort()
    end

    vars = @variables begin
        V_m(t) = V_m
        Phi(t) = Phi
    end

    equations = Equation[
        V_m ~ port_p.V_m - port_n.V_m,
        Phi ~ port_p.Phi,
        0 ~ port_p.Phi + port_n.Phi
    ]

    return System(equations, t, vars, pars; name, systems)
end
