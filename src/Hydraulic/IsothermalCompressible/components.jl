"""
    Source(;name, p)

Fixed pressure source

# Parameters:

  - `p`: [Pa] set pressure

# Connectors:

  - `port`: hydraulic port
"""
function Source(;name, p)
    pars = @parameters begin
        p = p
    end

    vars = []
    
    systems = @named begin
        port = HydraulicPort(; p_int = p)
    end
    
    eqs = [
        port.p ~ p
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end



function InputSource(;name, p_int)
    pars = @parameters begin
        p_int = p_int
    end

    vars = []
    
    systems = @named begin
        port = HydraulicPort(; p_int)
        input = RealInput()
    end
    
    eqs = [
        port.p ~ input.u
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end


"""
    FixedVolume(fluid; name, vol, p_int)

fixed fluid volume where `fluid` specifies the medium

# Parameters:

  - `vol`: [m^3] fixed volume
  - `p_int`: [Pa] Initial pressure
  

# Connectors:

  - `port`: hydraulic port
"""
function FixedVolume(; name, fluid, vol, p_int)
    pars = @parameters begin
        vol = vol
        p_int = p_int
        fluid=fluid
    end

    vars = @variables begin
        rho(t) = density(fluid, p_int)
        drho(t) = 0
    end
    
    systems = @named begin
        port = HydraulicPort(; p_int)
    end

    # let -------------------
    dm = port.dm
    
    eqs = [
        D(rho) ~ drho
        rho ~ density(fluid, port.p)
        
        dm ~ drho*vol 
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    LaminarResistance(fluid, shape=:circle; name, p_int, area, length, perimeter=2*sqrt(area*pi))

laminar pipe resistance

# Parameters:

  - `vol`: [m^3] fixed volume
  - `p_int`: [Pa] Initial pressure
  

# Connectors:

  - `port`: hydraulic port
"""
function PipeBase(;name, fluid, shape=Shapes.circle, p_int, area, length, perimeter=2*sqrt(area*pi))
    pars = @parameters begin
        p_int = p_int
        area = area
        length = length
        perimeter = perimeter
        fluid=fluid
        shape=shape
    end

    vars = []
    
    systems = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
    end
    

    # let ----------------------
    Δp = port_a.p - port_b.p
    dm = port_a.dm
    p = (port_a.p + port_b.p)/2

    d_h = 4*area/perimeter

    ρ = density(fluid, p)
    μ = viscosity(fluid)
    Φ = shape_factor(shape)
    f = friction_factor(dm, area, d_h, ρ, μ, Φ)
    u = dm/(ρ*area)

    eqs = [        
        Δp ~ 1/2 * ρ * u^2 * f * (length/d_h)
        0 ~ port_a.dm + port_b.dm
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

function Pipe(N; name, fluid, shape=Shapes.circle, p_int, area, length, perimeter=2*sqrt(area*pi))

    @assert(N>1, "the pipe component must be defined with more than 1 segment (i.e. N>1), found N=$N")

    pars = @parameters begin
        p_int = p_int
        area = area
        length = length
        perimeter = perimeter
        fluid=fluid
        shape=shape
    end
    
    vars = []
    
    ports = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
    end


    pipe_bases = []
    for i=1:N-1
        x = PipeBase(;name=Symbol("p$i"), fluid=ParentScope(fluid), shape=ParentScope(shape), p_int=ParentScope(p_int), area=ParentScope(area), length=ParentScope(length)/(N-1), perimeter=ParentScope(perimeter))
        push!(pipe_bases, x)
    end

    volumes = []
    for i=1:N
        x = FixedVolume(; name=Symbol("v$i"), fluid=ParentScope(fluid), vol=ParentScope(area)*ParentScope(length)/N, p_int=ParentScope(p_int))
        push!(volumes, x)
    end
    

    eqs = [
        connect(volumes[1].port, pipe_bases[1].port_a, port_a)
        connect(volumes[end].port, pipe_bases[end].port_b, port_b)
    ]

    for i=2:N-1
        eq = connect(volumes[i].port, pipe_bases[i-1].port_b, pipe_bases[i].port_a)
        push!(eqs, eq)
    end


    ODESystem(eqs, t, vars, pars; name, systems=[ports; pipe_bases; volumes])
end


function Actuator(direction; name, fluid, p_int, x_int, area, dead_volume)

    pars = @parameters begin
        fluid = fluid
        p_int = p_int
        x_int = x_int
        area = area
        dead_volume = dead_volume
    end

    vars = @variables begin
        x(t) = x_int
        dx(t) = 0
        rho(t) = density(fluid, p_int)
        drho(t) = 0
    end

    systems = @named begin
        port = HydraulicPort(; p_int)
        flange = MechanicalPort()
    end

    # let -------------
    vol = dead_volume + area*x

    eqs = [
        D(x) ~ dx
        D(rho) ~ drho

        dx ~ flange.v*direction
        rho ~ density(fluid, port.p)
        
        port.dm ~ drho*vol + rho*area*dx
        flange.f ~ -port.p*area*direction
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

