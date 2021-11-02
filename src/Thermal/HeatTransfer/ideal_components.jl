function ThermalGround(; name)
    @named hp = HeatPort()
    eqs = [hp.T ~ 0]
    ODESystem(eqs, t, systems=[hp], name=name)
end

function HeatCapacitor(; name, C=1.0)
    c_th = C
    
    @named hp = HeatPort()
    @parameters C
    @variables T(t) dt(t)

    D = Differential(t)
    eqs = [
        T ~ hp.T
        dt ~ D(T)
        D(T) ~ hp.Q_flow / C
        ]
    ODESystem(eqs, t, [T, dt], [C], systems=[hp], defaults=Dict(C => c_th), name=name)
end

function ThermalConductor(; name, G=1.0)
    g_th = G

    @named hp1 = HeatPort()
    @named hp2 = HeatPort()
    @parameters G
    @variables Q_flow(t) T(t)

    eqs = [
        T ~ hp1.T - hp2.T
        Q_flow ~ G*T
        Q_flow ~ hp1.Q_flow
        -Q_flow ~ hp2.Q_flow
    ]
    ODESystem(eqs, t, [Q_flow, T], [G], systems=[hp1, hp2], defaults=Dict(G => g_th), name=name)
end

function ThermalResistor(; name, R=1.0)
    r_th = R
    @named hp1 = HeatPort()
    @named hp2 = HeatPort()
    @parameters R
    @variables Q_flow(t) T(t)

    eqs = [
        T ~ R*Q_flow
        T ~ hp1.T - hp2.T
        hp1.Q_flow ~ Q_flow
        hp2.Q_flow ~ -Q_flow
    ]
    ODESystem(eqs, t, [Q_flow, T], [R], systems=[hp1, hp2], defaults=Dict(R => r_th), name=name)
end

function ConvectiveConductor(; name, G=1.0)
    g_c = G

    @named solidport = HeatPort()
    @named fluidport = HeatPort()
    @parameters G # Convective thermal conductance
    @variables Q_flow(t) dT(t) 

    eqs = [
        dT ~ solidport.T - fluidport.T
        solidport.Q_flow ~ Q_flow
        fluidport.Q_flow ~ -Q_flow
        dT ~ G*Q_flow
    ]
    ODESystem(eqs, t, [Q_flow, dT], [G], systems=[solidport, fluidport], defaults=Dict(G => g_c), name=name)
end

function ConvectiveResistor(; name, R=1.0)
    r_c = R

    @named solidport = HeatPort()
    @named fluidport = HeatPort()
    @parameters R # Convective thermal resistance
    @variables Q_flow(t) dT(t)

    eqs = [
        dT ~ solidport.T - fluidport.T
        solidport.Q_flow ~ Q_flow
        fluidport.Q_flow ~ -Q_flow
        dT ~ R*Q_flow
    ]
    ODESystem(eqs, t, [Q_flow, dT], [R], systems=[solidport, fluidport], defaults=Dict(R => r_c), name=name)
end

function BodyRadiation(; name, G=1.0)
    g_r = G 
    σ   = 5.6703744191844294e-8 # Stefan-Boltzmann constant
    
    @named hp1 = HeatPort()
    @named hp2 = HeatPort()
    @parameters G # Net radiation conductance between two surfaces 
    @variables Q_flow(t)

    eqs = [
        Q_flow ~ G*σ*(hp1.T^4 - hp2.T^4)
    ]
    ODESystem(eqs, t, [Q_flow], [G], systems=[hp1, hp2], defaults=Dict(G => g_r), name=name)
end

function ThermalCollector(; name, N=1)
    hp = []
    for i in 1:N
        _hp = HeatPort(name=Symbol(:hp, i))
        push!(hp, _hp)
    end
    @named collector_port = HeatPort()

    eqs = [
        collector_port.Q_flow + sum(k -> k.Q_flow, hp) ~ 0
        collector_port.T ~ hp[1].T
    ]
    for i in 1:N-1
        push!(eqs, hp[i].T ~ hp[i+1].T)
    end
    ODESystem(eqs, t, [], [], systems=[hp..., collector_port], name=name)
end
