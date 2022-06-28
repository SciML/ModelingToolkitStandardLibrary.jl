"""
    HeatCapacitor(; name, C, T_start=273.15 + 20)    

Lumped thermal element storing heat

# States:
- `T`: [`K`] Temperature of element
- `der_T`: [`K/s`] Time derivative of temperature

# Connectors:
- `port`

# Parameters:
- `C`: [`J/K`] Heat capacity of element (= cp*m)
- `T_start`: [`K`] Initial temperature of element
"""
function HeatCapacitor(; name, C, T_start = 273.15 + 20)
    @named port = HeatPort()
    @parameters C = C
    sts = @variables begin
        T(t) = T_start
        der_T(t) = 0.0
    end

    D = Differential(t)
    eqs = [T ~ port.T
           der_T ~ port.Q_flow / C
           D(T) ~ der_T]
    ODESystem(eqs, t, sts, [C]; systems = [port], name = name)
end

"""
    ThermalConductor(;name, G) 

Lumped thermal element transporting heat without storing it.

# States:
see [`Element1D`](@ref)

# Connectors:
`port_a`
`port_b`

# Parameters:
- `G`: [`W/K`] Constant thermal conductance of material
"""
function ThermalConductor(; name, G)
    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d
    pars = @parameters G = G
    eqs = [
        Q_flow ~ G * dT,
    ]
    extend(ODESystem(eqs, t, [], pars; name = name), element1d)
end

"""
    ThermalResistor(; name, R) 

Lumped thermal element transporting heat without storing it.

# States:
- `dT`:  [`K`] Temperature difference across the component a.T - b.T
- `Q_flow`: [`W`] Heat flow rate from port a -> port b

# Connectors:
- `port_a`
- `port_b`

# Parameters:
- `R`: [`K/W`] Constant thermal resistance of material
"""
function ThermalResistor(; name, R)
    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d
    pars = @parameters R = R
    eqs = [
        dT ~ R * Q_flow,
    ]

    extend(ODESystem(eqs, t, [], pars; name = name), element1d)
end

"""
    ConvectiveConductor(; name, G)

Lumped thermal element for heat convection.

# States:
- `dT`:  [`K`] Temperature difference across the component solid.T - fluid.T
- `Q_flow`: [`W`] Heat flow rate from `solid` -> `fluid`

# Connectors:
- `solid`
- `fluid`

# Parameters:
- `G`: [W/K] Convective thermal conductance
"""
function ConvectiveConductor(; name, G)
    @named solid = HeatPort()
    @named fluid = HeatPort()
    @parameters G = G
    sts = @variables Q_flow(t)=0.0 dT(t)=0.0
    eqs = [dT ~ solid.T - fluid.T
           solid.Q_flow ~ Q_flow
           fluid.Q_flow ~ -Q_flow
           dT ~ G * Q_flow]
    ODESystem(eqs, t, sts, [G]; systems = [solid, fluid], name = name)
end

"""
    ConvectiveResistor(; name, R)

Lumped thermal element for heat convection.

# States:
- `dT`:  [`K`] Temperature difference across the component solid.T - fluid.T
- `Q_flow`: [`W`] Heat flow rate from `solid` -> `fluid`

# Connectors:
- `solid`
- `fluid`

# Parameters:
- `R`: [`K/W`] Constant thermal resistance of material
"""
function ConvectiveResistor(; name, R)
    @named solid = HeatPort()
    @named fluid = HeatPort()
    @parameters R = R
    sts = @variables Q_flow(t)=0.0 dT(t)=0.0
    eqs = [dT ~ solid.T - fluid.T
           solid.Q_flow ~ Q_flow
           fluid.Q_flow ~ -Q_flow
           dT ~ R * Q_flow]
    ODESystem(eqs, t, sts, [R]; systems = [solid, fluid], name = name)
end

"""
    BodyRadiation(; name, G)

Lumped thermal element for radiation heat transfer.

# States:
- `dT`:  [`K`] Temperature difference across the component a.T - b.T
- `Q_flow`: [`W`] Heat flow rate from port a -> port b

# Connectors:
- `port_a`
- `port_b`

# Parameters:
- `G`: [m^2] Net radiation conductance between two surfaces
"""
function BodyRadiation(; name, G)
    sigma = 5.6703744191844294e-8 # Stefan-Boltzmann constant TODO: extract into physical constants module or use existing one

    @named element1d = Element1D()
    @unpack Q_flow, dT = element1d
    @unpack port_a, port_b = element1d
    pars = @parameters G = G
    eqs = [
        Q_flow ~ G * sigma * (port_a.T^4 - port_b.T^4),
    ]

    extend(ODESystem(eqs, t, [], pars; name = name), element1d)
end

"""
    ThermalCollector(; name, m=1)

Collects `m` heat flows

This is a model to collect the heat flows from `m` heatports to one single heatport.

# States:

# Connectors:
- `port_a1` to `port_am`
- `port_b`

# Parameters:
- `m`: Number of heat ports (e.g. m=2: `port_a1`, `port_a2`)
"""
function ThermalCollector(; name, m::Integer = 1)
    port_a = [HeatPort(name = Symbol(:port_a, i)) for i in 1:m]
    @named port_b = HeatPort()
    eqs = [port_b.Q_flow + sum(k -> k.Q_flow, port_a) ~ 0
           port_b.T ~ port_a[1].T]
    for i in 1:(m - 1)
        push!(eqs, port_a[i].T ~ port_a[i + 1].T)
    end
    ODESystem(eqs, t, [], []; systems = [port_a..., port_b], name = name)
end
