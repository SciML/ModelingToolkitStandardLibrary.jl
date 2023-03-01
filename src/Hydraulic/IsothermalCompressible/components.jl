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

fixed fluid volume

# Parameters:

  - `vol`: [m^3] fixed volume
  - `p_int`: [Pa] Initial pressure
  

# Connectors:

  - `port`: hydraulic port
"""
function LaminarResistance(;fluid, shape=Int(circle), name, p_int, area, length, perimeter=2*sqrt(area*pi))
    pars = @parameters begin
        p_int = p_int
        area = area
        length = length
        perimeter = perimeter
        fluid=fluid
        shape=shape
    end

    vars = @variables begin
        v(t) = 0
    end
    
    systems = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
    end
    

    # let ----------------------
    Δp = port_a.p - port_b.p
    dm = port_a.dm

    d_h = 4*area/perimeter
    rho = density(fluid, port_a.p)
    Re = rho*v*d_h/viscosity(fluid)
    f = friction_factor(shape)/Re
    

    eqs = [        
        Δp ~ 1/2 * rho * v^2 * f * length
        dm ~ rho*v*area
        0 ~ port_a.dm - port_b.dm
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

