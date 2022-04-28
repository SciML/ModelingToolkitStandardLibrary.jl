"""
    HeatCapacitor(; name, C=1.0, T_start=293.15 + 20)    

Lumped thermal element storing heat

# Parameters:
- `C`: [J/K] Heat capacity of element (= cp*m)
- `T_start`: Initial temperature of element

# States:
- `T`: [K] Temperature of element
- `der_T`: [K/s] Time derivative of temperature
"""
function HeatCapacitor(; name, C=1.0, T_start=293.15 + 20)    
    @named port = HeatPort()
    @parameters C=C
    sts = @variables begin
        T(t)=T_start
        der_T(t)
    end

    D = Differential(t)
    eqs = [
        T ~ port.T
        der_T ~ D(T)
        D(T) ~ port.Q_flow / C
    ]
    ODESystem(eqs, t, sts, [C]; systems=[port], name=name)
end

"""
    ThermalConductor(;name, G=1.0) 

Lumped thermal element transporting heat without storing it.

# Parameters:
- `G`: [W/K] Constant thermal conductance of material
"""
function ThermalConductor(;name, G=1.0)   
    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d
    pars = @parameters G=G
    eqs = [
        Q_flow ~ G * dT
    ]
    extend(ODESystem(eqs, t, [], pars; name=name), element1d)
end

"""
    ThermalResistor(; name, R=1.0) 

Lumped thermal element transporting heat without storing it.

# Parameters:
- `R`: [K/W] Constant thermal resistance of material
"""
function ThermalResistor(; name, R=1.0)   
    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d
    pars = @parameters R=R
    eqs = [
        dT ~ R * Q_flow
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), element1d)
end

"""
    ConvectiveConductor(; name, G=1.0)

Lumped thermal element for heat convection.

# Parameters:
- `G`: [W/K] Convective thermal conductance

# States:
- `dT`:  [K] Temperature difference across the component solid.T - fluid.T
- `Q_flow`: [W] Heat flow rate from solid -> fluid
"""
function ConvectiveConductor(; name, G=1.0)
    @named solid = HeatPort()
    @named fluid = HeatPort()
    @parameters G=G
    sts = @variables Q_flow(t) dT(t)
    eqs = [
        dT ~ solid.T - fluid.T
        solid.Q_flow ~ Q_flow
        fluid.Q_flow ~ -Q_flow
        dT ~ G*Q_flow
    ]
    ODESystem(eqs, t, sts, [G]; systems=[solid, fluid], name=name)
end

"""
    ConvectiveResistor(; name, R=1.0)

Lumped thermal element for heat convection.

# Parameters:
- `R`: [K/W] Constant thermal resistance of material

# States:
- `dT`:  [K] Temperature difference across the component solid.T - fluid.T
- `Q_flow`: [W] Heat flow rate from solid -> fluid
"""
function ConvectiveResistor(; name, R=1.0)
    @named solidport = HeatPort()
    @named fluidport = HeatPort()
    @parameters R=R
    sts = @variables Q_flow(t) dT(t) 
    eqs = [
        dT ~ solidport.T - fluidport.T
        solidport.Q_flow ~ Q_flow
        fluidport.Q_flow ~ -Q_flow
        dT ~ R*Q_flow
    ]
    ODESystem(eqs, t, sts, [R]; systems=[solidport, fluidport], name=name)
end

"""
    BodyRadiation(; name, G=1.0)

Lumped thermal element for radiation heat transfer.

# Parameters:
- `G`: [m^2] Net radiation conductance between two surfaces
"""
function BodyRadiation(; name, G=1.0)
    sigma = 5.6703744191844294e-8 # Stefan-Boltzmann constant TODO: extract into physical constants module or use existing one

    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d
    @unpack port_a, port_b = element1d
    pars = @parameters G=G
    eqs = [
        Q_flow ~ G * sigma * (port_a.T^4 - port_b.T^4)
    ]
    
    extend(ODESystem(eqs, t, [], pars; name=name), element1d)
end

"""
    ThermalCollector(; name, m=1)

Collects m heat flows

This is a model to collect the heat flows from `m` heatports to one single heatport.
# Parameters:
- `m`: Number of heat ports (e.g. m=2: `port_a1`, `port_a2`)
"""
function ThermalCollector(; name, m=1)
    port_a = [HeatPort(name=Symbol(:port_a, i)) for i in 1:m]
    @named port_b = HeatPort()
    eqs = [
        port_b.Q_flow + sum(k -> k.Q_flow, port_a) ~ 0
        port_b.T ~ port_a[1].T
    ]
    for i in 1:m-1
        push!(eqs, port_a[i].T ~ port_a[i+1].T)
    end
    ODESystem(eqs, t, [], []; systems=[port_a..., port_b], name=name)
end
