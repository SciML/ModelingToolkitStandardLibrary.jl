"""
    Ground(; name)

Zero magnetic potential.
"""
@component function Ground(; name)
    pars = @parameters begin
    end

    systems = @named begin
        port = PositiveMagneticPort()
    end

    vars = @variables begin
    end

    equations = Equation[
        port.V_m ~ 0,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    Idle(;name)

Idle running branch.
"""
@component function Idle(; name)
    @named two_port = TwoPort()
    @unpack Phi = two_port

    pars = @parameters begin
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        Phi ~ 0,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, two_port)
end

"""
    Short(;name)

Short cut branch.
"""
@component function Short(; name)
    @named two_port = TwoPort()
    @unpack V_m = two_port

    pars = @parameters begin
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        V_m ~ 0,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, two_port)
end

"""
    Crossing(;name)

Crossing of two branches.

This is a simple crossing of two branches. The ports port_p1 and port_p2 are connected, as well as port_n1 and port_n2.
"""
@component function Crossing(; name)
    pars = @parameters begin
    end

    systems = @named begin
        port_p1 = PositiveMagneticPort()
        port_p2 = PositiveMagneticPort()
        port_n1 = NegativeMagneticPort()
        port_n2 = NegativeMagneticPort()
    end

    vars = @variables begin
    end

    equations = Equation[
        connect(port_p1, port_p2),
        connect(port_n1, port_n2),
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    ConstantPermeance(; name, G_m = 1.0)

Constant permeance.

# Parameters:

  - `G_m`: [H] Magnetic permeance
"""
@component function ConstantPermeance(; name, G_m = 1.0)
    @named two_port = TwoPort()
    @unpack V_m, Phi = two_port

    pars = @parameters begin
        G_m = G_m, [description = "Magnetic permeance"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        Phi ~ G_m * V_m,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, two_port)
end

"""
    ConstantReluctance(; name, R_m = 1.0)

Constant reluctance.

# Parameters:

  - `R_m`: [H^-1] Magnetic reluctance
"""
@component function ConstantReluctance(; name, R_m = 1.0)
    @named two_port = TwoPort(; Phi = 0.0)
    @unpack V_m, Phi = two_port

    pars = @parameters begin
        R_m = R_m, [description = "Magnetic reluctance"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        V_m ~ Phi * R_m,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, two_port)
end

"""
    ElectroMagneticConverter(; name, N, Phi)

Ideal electromagnetic energy conversion.

The electromagnetic energy conversion is given by Ampere's law and Faraday's law respectively
V_m = N * i
N * dΦ/dt = -v

Initial magnetic flux flowing into the port_p can be set with `Phi` ([Wb])

# Parameters:

  - `N`: Number of turns
"""
@component function ElectroMagneticConverter(; name, N = nothing, Phi = nothing, v = nothing, i = nothing)
    @named two_port = TwoPort(; Phi)
    @unpack V_m, Phi = two_port

    pars = @parameters begin
        N = N, [description = "Number of turns"]
    end

    systems = @named begin
        p = Pin()
        n = Pin()
    end

    vars = @variables begin
        v(t) = v
        i(t) = i
    end

    equations = Equation[
        v ~ p.v - n.v,
        0 ~ p.i + n.i,
        i ~ p.i,
        #converter equations:
        V_m ~ i * N, # Ampere's law
        D(Phi) ~ -v / N,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, two_port)
end

"""
    EddyCurrent(;name, Phi, rho = 0.098e-6, l = 1, A = 1)

For modelling of eddy current in a conductive magnetic flux tube.
Initial magnetic flux flowing into the port_p can be set with `Phi` ([`Wb`])

# Parameters:

  - `rho`: [ohm * m] Resistivity of flux tube material (default: Iron at 20degC)
  - `l`: [m] Average length of eddy current path
  - `A`: [m^2] Cross sectional area of eddy current path
"""
@component function EddyCurrent(; name, Phi = nothing, rho = 0.098e-6, l = 1, A = 1)
    @named two_port = TwoPort(; Phi)
    @unpack V_m, Phi = two_port

    R = rho * l / A # Electrical resistance of eddy current path

    pars = @parameters begin
        rho = rho, [description = "Resistivity of flux tube material"]
        l = l, [description = "Average length of eddy current path"]
        A = A, [description = "Cross sectional area of eddy current path"]
        R = R
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = Equation[
        D(Phi) ~ V_m * R,
    ]

    sys = System(equations, t, vars, pars; name, systems)
    return extend(sys, two_port)
end
