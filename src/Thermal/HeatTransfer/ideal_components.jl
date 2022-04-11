function ThermalGround(; name)
    @named a = HeatPort()
    eqs = [a.T ~ 0]
    ODESystem(eqs, t, systems=[a], name=name)
end

"""
Lumped thermal element storing heat
"""
function HeatCapacitor(; name, 
    C=1.0, # [J/K] Heat capacity of element
    )    
    @named a = HeatPort()
    @parameters C=C
    sts = @variables begin
        T(t) # Temperature of element
        der_T(t) # "Time derivative of temperature
    end

    D = Differential(t)
    eqs = [
        T ~ a.T
        der_T ~ D(T)
        D(T) ~ a.Q_flow / C
        ]
    ODESystem(eqs, t, sts, [C]; systems=[a], name=name)
end

"""
Lumped thermal element transporting heat without storing it.
"""
function ThermalConductor(;name, 
    G=1.0, # [W/K] Constant thermal conductance of material
    )   
    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d
    pars = @parameters G=G
    eqs = [
        Q_flow ~ G * dT
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), element1d)
end

"""
Lumped thermal element transporting heat without storing it.
"""
function ThermalResistor(; name,
    R=1.0, # [K/W] Constant thermal resistance of material
    )   
    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d
    pars = @parameters R=R
    eqs = [
        dT ~ R * Q_flow
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), element1d)
end

"""
Lumped thermal element for heat convection.
"""
function ConvectiveConductor(; name, 
    G=1.0, # [W/K] Convective thermal conductance
    )
    @named solidport = HeatPort()
    @named fluidport = HeatPort()
    @parameters G=G
    sts = @variables begin
        Q_flow(t) # [W] Heat flow rate from solid -> fluid
        dT(t) # [K] Temperature difference solid.T - fluid.T
    end

    eqs = [
        dT ~ solidport.T - fluidport.T
        solidport.Q_flow ~ Q_flow
        fluidport.Q_flow ~ -Q_flow
        dT ~ G*Q_flow
    ]
    ODESystem(eqs, t, sts, [G]; systems=[solidport, fluidport], name=name)
end

"""
Lumped thermal element for heat convection.
"""
function ConvectiveResistor(; name, 
    R=1.0, # [K/W] Convective thermal resistance
    )
    @named solidport = HeatPort()
    @named fluidport = HeatPort()
    @parameters R=R
    sts = @variables begin
        Q_flow(t) # [W] Heat flow rate from solid -> fluid
        dT(t) # [K] Temperature difference solid.T - fluid.T
    end

    eqs = [
        dT ~ solidport.T - fluidport.T
        solidport.Q_flow ~ Q_flow
        fluidport.Q_flow ~ -Q_flow
        dT ~ R*Q_flow
    ]
    ODESystem(eqs, t, sts, [R]; systems=[solidport, fluidport], name=name)
end

"""
Lumped thermal element for radiation heat transfer.
"""
function BodyRadiation(; name, 
    G=1.0, # [m^2] Net radiation conductance between two surfaces
    )
    sigma = 5.6703744191844294e-8 # Stefan-Boltzmann constant

    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d
    pars = @parameters G=G
    eqs = [
        Q_flow ~ G * sigma * (element1d.a.T^4 - element1d.b.T^4)
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), element1d)
end

"""
This is a model to collect the heat flows from `N` heatports to one single heatport.
"""
function ThermalCollector(; name, N=1)
    hp = [HeatPort(name=Symbol(:hp, i)) for i in 1:N]
    @named collector_port = HeatPort()
    eqs = [
        collector_port.Q_flow + sum(k -> k.Q_flow, hp) ~ 0
        collector_port.T ~ hp[1].T
    ]
    for i in 1:N-1
        push!(eqs, hp[i].T ~ hp[i+1].T)
    end
    ODESystem(eqs, t, [], []; systems=[hp..., collector_port], name=name)
end
