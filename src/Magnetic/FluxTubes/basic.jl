"""
    Ground(;name)

Zero magnetic potential.
"""
@component function Ground(; name)
    @named port = PositiveMagneticPort()
    eqs = [port.V_m ~ 0]
    ODESystem(eqs, t, [], [], systems = [port], name = name)
end

"""
    Idle(;name)

Idle running branch.
"""
@component function Idle(; name)
    @named two_port = TwoPort()
    @unpack Phi = two_port
    eqs = [
        Phi ~ 0,
    ]
    extend(ODESystem(eqs, t, [], [], systems = [], name = name), two_port)
end

"""
    Short(;name)

Short cut branch.
"""
@component function Short(; name)
    @named two_port = TwoPort()
    @unpack V_m = two_port
    eqs = [
        V_m ~ 0,
    ]
    extend(ODESystem(eqs, t, [], [], systems = [], name = name), two_port)
end

"""
    Crossing(;name)

Crossing of two branches.

This is a simple crossing of two branches. The ports port_p1 and port_p2 are connected, as well as port_n1 and port_n2.
"""
@component function Crossing(; name)
    @named port_p1 = PositiveMagneticPort()
    @named port_p2 = PositiveMagneticPort()
    @named port_n1 = NegativeMagneticPort()
    @named port_n2 = NegativeMagneticPort()
    eqs = [
        connect(port_p1, port_p2),
        connect(port_n1, port_n2),
    ]
    ODESystem(eqs, t, [], [], systems = [port_p1, port_p2, port_n1, port_n2], name = name)
end

"""
    ConstantPermeance(;name, G_m=1.0)

Constant permeance.

# Parameters:

  - `G_m`: [H] Magnetic permeance
"""
@component function ConstantPermeance(; name, G_m = 1.0)
    @named two_port = TwoPort()
    @unpack V_m, Phi = two_port
    @parameters G_m = G_m
    eqs = [
        Phi ~ G_m * V_m,
    ]
    extend(ODESystem(eqs, t, [], [G_m], name = name), two_port)
end

"""
    ConstantReluctance(;name, R_m=1.0)

Constant reluctance.

# Parameters:

  - `R_m`: [H^-1] Magnetic reluctance
"""
@component function ConstantReluctance(; name, R_m = 1.0)
    @named two_port = TwoPort()
    @unpack V_m, Phi = two_port
    @parameters R_m = R_m
    eqs = [
        V_m ~ Phi * R_m,
    ]
    extend(ODESystem(eqs, t, [], [R_m], name = name), two_port)
end

"""
    ElectroMagneticConverter(;name, N, Phi_start=0.0)

Ideal electromagnetic energy conversion.

The electromagnetic energy conversion is given by Ampere's law and Faraday's law respectively
V_m = N * i
N * dΦ/dt = -v

# Parameters:

  - `N`: Number of turns
  - `Phi_start`: [Wb] Initial magnetic flux flowing into the port_p
"""
@component function ElectroMagneticConverter(; name, N, Phi_start = 0.0)
    @named port_p = PositiveMagneticPort()
    @named port_n = NegativeMagneticPort()
    @named p = Pin()
    @named n = Pin()

    sts = @variables v(t) i(t) V_m(t) Phi(t)=Phi_start
    pars = @parameters N = N
    eqs = [v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
        V_m ~ port_p.V_m - port_n.V_m
        0 ~ port_p.Phi + port_n.Phi
        Phi ~ port_p.Phi
    #converter equations:
        V_m ~ i * N # Ampere's law
        D(Phi) ~ -v / N]
    ODESystem(eqs, t, sts, pars, systems = [port_p, port_n, p, n], name = name)
end

"""
    EddyCurrent(;name, rho=0.098e-6, l=1, A=1, Phi_start=0.0)

For modelling of eddy current in a conductive magnetic flux tube.

# Parameters:

  - `rho`: [ohm * m] Resistivity of flux tube material (default: Iron at 20degC)
  - `l`: [m] Average length of eddy current path
  - `A`: [m^2] Cross sectional area of eddy current path
  - `Phi_start`: [Wb] Initial magnetic flux flowing into the port_p
"""
@component function EddyCurrent(; name, rho = 0.098e-6, l = 1, A = 1, Phi_start = 0.0)
    @named two_port = TwoPort(Phi_start = Phi_start)
    @unpack V_m, Phi = two_port
    @parameters R = rho * l / A # Electrical resistance of eddy current path
    eqs = [
        D(Phi) ~ V_m * R,
    ]
    extend(ODESystem(eqs, t, [], [R], name = name), two_port)
end
