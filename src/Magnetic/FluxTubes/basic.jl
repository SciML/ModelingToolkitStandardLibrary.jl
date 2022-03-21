function Ground(;name)
    @named port = PositiveMagneticPort()
    eqs = [port.V_m ~ 0]
    ODESystem(eqs, t, [], [], systems=[port], name=name)
end

function Idle(;name)
    @named two_port = TwoPort()
    @unpack Phi = two_port
    eqs = [
        Phi ~ 0,
    ]
    extend(ODESystem(eqs, t, [Phi], [], systems=[], name=name), two_port)
end

function Short(;name)
    port_p = PositiveMagneticPort()
    port_n = NegativeMagneticPort()
    eqs = [
        connect(port_p, port_n),
    ]
    ODESystem(eqs, t, [], [], systems=[port_p, port_n], name=name)
end

function Crossing(;name)
    port_p1 = PositiveMagneticPort()
    port_p2 = PositiveMagneticPort()
    port_n1 = NegativeMagneticPort()
    port_n2 = NegativeMagneticPort()
    eqs = [
        connect(port_p1, port_p2),
        connect(port_n1, port_n2),
    ]
    ODESystem(eqs, t, [], [], systems=[port_p1, port_p2, port_n1, port_n2], name=name)
end

function ConstantPermeance(;name, G_m=1.0)
    val = G_m
    @named two_port = TwoPort()
    @unpack V_m, Phi = two_port
    @variables G_m(t)
    eqs = [
        Phi ~ G_m * V_m,
    ]
    extend(ODESystem(eqs, t, [V_m, Phi, G_m], [], systems=[], defaults=Dict(G_m => val), name=name), two_port)
end

function ConstantReluctance(;name, R_m=1.0)
    val = R_m
    @named two_port = TwoPort()
    @unpack V_m, Phi = two_port
    @variables R_m(t)
    eqs = [
        V_m ~ Phi * R_m,
    ]
    extend(ODESystem(eqs, t, [V_m, Phi, R_m], [], systems=[], defaults=Dict(R_m => val), name=name), two_port)
end
