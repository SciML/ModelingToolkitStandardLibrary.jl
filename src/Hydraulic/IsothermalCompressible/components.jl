"""
    Source(; p, name)

Fixed pressure source

# Parameters:
- `p`: [Pa] set pressure (set by `p` argument)

# Connectors:
- `port`: hydraulic port
"""
@component function Source(; p, name)
    pars = @parameters begin p = p end

    vars = []

    systems = @named begin port = HydraulicPort(; p_int = p) end

    eqs = [
        port.p ~ p,
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    InputSource(; p_int, name)

Fixed pressure source

# Parameters:
- `p_int`: [Pa] initial pressure (set by `p_int` argument)

# Connectors:
- `port`: hydraulic port
"""
@component function InputSource(; p_int, name)
    pars = @parameters begin p_int = p_int end

    vars = []

    systems = @named begin
        port = HydraulicPort(; p_int)
        input = RealInput()
    end

    eqs = [
        port.p ~ input.u,
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    Cap(; p_int, name)

Caps a hydrualic port to prevent mass flow in or out.

# Parameters:
- `p_int`: [Pa] initial pressure (set by `p_int` argument)

# Connectors:
- `port`: hydraulic port
"""
@component function Cap(; p_int, name)
    pars = @parameters p_int = p_int

    vars = @variables p(t) = p_int

    systems = @named begin port = HydraulicPort(; p_int = p_int) end

    eqs = [port.p ~ p
           port.dm ~ 0]

    ODESystem(eqs, t, vars, pars; name, systems)
end


@component function Open(; p_int, name)
    pars = @parameters p_int = p_int

    vars = @variables begin
        p(t) = p_int
        dm(t) = 0
    end
        

    systems = @named begin port = HydraulicPort(; p_int = p_int) end

    eqs = [port.p ~ p
           port.dm ~ dm]

    ODESystem(eqs, t, vars, pars; name, systems)
end


"""
    FixedVolume(; vol, p_int, name)

Fixed fluid volume.

# Parameters:
- `vol`: [m^3] fixed volume
- `p_int`: [Pa] initial pressure

# Connectors:
- `port`: hydraulic port
"""
@component function FixedVolume(; vol, p_int, name)
    pars = @parameters begin
        vol = vol
        p_int = p_int
    end

    systems = @named begin port = HydraulicPort(; p_int) end

    vars = @variables begin
        rho(t) = liquid_density(port)
        drho(t) = 0
    end

    # let -------------------
    dm = port.dm

    eqs = [D(rho) ~ drho
           rho ~ full_density(port)
           dm ~ drho * vol]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    TubeBase(; p_int, area, length, perimeter=2*sqrt(area*pi), shape_factor=64, name)

Internal flow model of the fully developed flow friction, ignoring any compressibility.

# Parameters:
- `p_int`: [Pa] initial pressure
- `area`: [m^2] tube cross sectional area
- `length`: [m] length of the pipe
- `perimeter`: [m] perimeter of the pipe cross section (needed only for non-circular pipes)
- `Φ`: shape factor, see `friction_factor` function (set by optional `shape_factor` argument, needed only for non-circular pipes).  

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
"""
@component function TubeBase(; p_int, area, length, perimeter = 2 * sqrt(area * pi),
                             shape_factor = 64, name)
    pars = @parameters begin
        p_int = p_int
        area = area
        length = length
        perimeter = perimeter
        Φ = shape_factor
    end

    vars = []

    systems = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
    end

    # let ----------------------
    Δp = port_a.p - port_b.p
    dm = port_a.dm

    d_h = 4 * area / perimeter

    ρ = (full_density(port_a) + full_density(port_b)) / 2
    μ = viscosity(port_a)

    f = friction_factor(dm, area, d_h, ρ, μ, Φ)
    u = dm / (ρ * area)

    eqs = [Δp ~ 1 / 2 * ρ * u^2 * f * (length / d_h)
           0 ~ port_a.dm + port_b.dm]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    Tube(N; p_int, area, length, effective_length=length, perimeter = 2 * sqrt(area * pi), shape_factor = 64, name)

Tube modeled with `N` segements which models the fully developed flow friction and compressibility.

# Parameters:
- `p_int`: [Pa] initial pressure 
- `area`: [m^2] tube cross sectional area 
- `length`: [m] real length of the tube 
- `effective_length`: [m] tube length to account for flow development and other restrictions, used to set the `TubeBase` length which calculates the flow resistance
- `perimeter`: [m] perimeter of the tube cross section (set by optional `perimeter` argument, needed only for non-circular tubes)
- `Φ`: shape factor, see `friction_factor` function (set by optional `shape_factor` argument, needed only for non-circular tubes)

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
"""
@component function Tube(N; p_int, area, length, effective_length = length,
                         perimeter = 2 * sqrt(area * pi),
                         shape_factor = 64, name)
    @assert(N>1,
            "the Tube component must be defined with more than 1 segment (i.e. N>1), found N=$N")

    #TODO: How to set an assert effective_length >= length ??
    pars = @parameters begin
        p_int = p_int
        area = area
        length = length
        effective_length = effective_length
        perimeter = perimeter
        Φ = shape_factor
    end

    vars = []

    ports = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
    end

    pipe_bases = []
    for i in 1:(N - 1)
        x = TubeBase(; name = Symbol("p$i"), shape_factor = ParentScope(Φ),
                     p_int = ParentScope(p_int), area = ParentScope(area),
                     length = ParentScope(effective_length) / (N - 1),
                     perimeter = ParentScope(perimeter))
        push!(pipe_bases, x)
    end

    volumes = []
    for i in 1:N
        x = FixedVolume(; name = Symbol("v$i"),
                        vol = ParentScope(area) * ParentScope(length) / N,
                        p_int = ParentScope(p_int))
        push!(volumes, x)
    end

    eqs = [connect(volumes[1].port, pipe_bases[1].port_a, port_a)
           connect(volumes[end].port, pipe_bases[end].port_b, port_b)]

    for i in 2:(N - 1)
        eq = connect(volumes[i].port, pipe_bases[i - 1].port_b, pipe_bases[i].port_a)
        push!(eqs, eq)
    end

    ODESystem(eqs, t, vars, pars; name, systems = [ports; pipe_bases; volumes])
end

"""
    FlowDivider(;p_int, n, name)

Reduces the flow from `port_a` to `port_b` by `n`.  Useful for modeling parallel tubes efficiently by placing a `FlowDivider` on each end of a tube.

# Parameters:
- `p_int`: [Pa] initial pressure 
- `n`: divide flow from `port_a` to `port_b` by `n`

# Connectors:
- `port_a`: full flow hydraulic port
- `port_b`: part flow hydraulic port
"""
@component function FlowDivider(; p_int, n, name)

    #TODO: assert n >= 1

    pars = @parameters begin
        n = n
        p_int = p_int
    end

    vars = @variables begin
        dm_a(t) = 0
        dm_b(t) = 0
    end

    systems = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
        open = Open(; p_int)

    end

    eqs = [
        connect(port_a, port_b, open.port)
        dm_a ~ port_a.dm
        dm_b ~ dm_a/n
        open.dm ~ dm_a - dm_b # extra flow dumps into an open port
        # port_b.dm ~ dm_b # divided flow goes to port_b
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    DynamicVolume(; p_int, x_int=0, area, dead_volume=0, direction=+1, minimum_volume=0, name)

Volume with moving wall.  The `direction` argument aligns the mechanical port with the hydraulic port, useful when connecting two dynamic volumes together in oppsing directions to create an actuator.
```
     ┌─────────────────┐ ───
     │                 │  ▲
                       │  │
dm ────►  dead volume  │  │ area
                       │  │  
     │                 │  ▼
     └─────────────────┤ ───
                       │
                       └─► x (= flange.v * direction)
```
# States:
- `x(t)`: [m] moving wall position
- `dx(t)`: [m/s] moving wall velocity
- `rho(t)`: [kg/m^3] density
- `drho(t)`: [kg/s-m^3] density derivative
- `vol(t)`: [m^3] volume
- `p(t)`: [Pa] dynamic pressure

# Parameters:
- `p_int`: [Pa] initial pressure
- `x_int`: [m] initial position of the moving wall
- `area`: [m^2] moving wall area
- `dead_volume`: [m^3] perimeter of the pipe cross section
- `minimum_volume`: [m^3] if `x`*`area` <= `minimum_volume` then mass flow `dm` shuts off

# Connectors:
- `port`: hydraulic port
- `flange`: mechanical translational port
"""
@component function DynamicVolume(; p_int, x_int = 0, area, dead_volume = 0, direction = +1,
                                  minimum_volume = 0, name)
    @assert (direction == +1)||(direction == -1) "direction arument must be +/-1, found $direction"

    pars = @parameters begin
        p_int = p_int
        x_int = x_int
        area = area
        dead_volume = dead_volume
    end

    systems = @named begin
        port = HydraulicPort(; p_int)
        flange = MechanicalPort()
    end

    vars = @variables begin
        x(t) = x_int
        dx(t) = 0
        rho(t) = liquid_density(port)
        drho(t) = 0
        vol(t) = dead_volume + area * x_int
        p(t) = p_int # p represents the dynamic pressure
    end

    # let
    dm = port.dm
    u = dm / (rho * area)

    eqs = [0 ~ IfElse.ifelse(vol >= minimum_volume, (p) - (port.p - (1 / 2) * rho * u^2),
                             dm)
           vol ~ dead_volume + area * x
           D(x) ~ dx
           D(rho) ~ drho
           dx ~ flange.v * direction
           rho ~ full_density(port)
           dm ~ drho * vol + rho * area * dx
           flange.f ~ -p * area * direction]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [flange.v => 0])
end

"""
    Valve(; p_a_int, p_b_int, area_int, Cd, name)

Valve with area input and discharge coefficient `Cd` defined by https://en.wikipedia.org/wiki/Discharge_coefficient

# Parameters:
- `p_a_int`: [Pa] initial pressure for `port_a`
- `p_b_int`: [Pa] initial pressure for `port_b`
- `area_int`: [m^2] initial valve opening
- `Cd`: discharge coefficient 

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
- `input`: real input setting the valve `area`.  Note: absolute value taken
"""
@component function Valve(reversible = false; p_a_int, p_b_int, area_int, Cd, name)
    pars = @parameters begin
        p_a_int = p_a_int
        p_b_int = p_b_int
        area_int = area_int
        Cd = Cd
    end

    systems = @named begin
        port_a = HydraulicPort(; p_int = p_a_int)
        port_b = HydraulicPort(; p_int = p_b_int)
        input = RealInput()
    end

    vars = []

    # let
    ρ = (full_density(port_a) + full_density(port_b)) / 2
    Δp = port_a.p - port_b.p
    dm = port_a.dm
    area = if reversible
        input.u
    else
        ifelse(input.u > 0, input.u , 0) 
    end

    eqs = [sign(Δp) * dm ~ sqrt(2 * abs(Δp) * ρ / Cd) * area
           0 ~ port_a.dm + port_b.dm]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [input.u => area_int])
end

@component function SpoolValve(reversible = false; p_a_int, p_b_int, x_int, Cd, d, name)
    pars = @parameters begin
        p_a_int = p_a_int
        p_b_int = p_b_int
        d = d
        x_int = x_int
        Cd = Cd
    end

    systems = @named begin
        port_a = HydraulicPort(; p_int = p_a_int)
        port_b = HydraulicPort(; p_int = p_b_int)
        flange = MechanicalPort()
        valve = Valve(reversible; p_a_int, p_b_int,
                      area_int = ParentScope(x_int) * 2π * ParentScope(d), Cd)
    end

    vars = @variables begin
        x(t) = x_int
        dx(t) = 0
    end

    eqs = [D(x) ~ dx
           flange.v ~ dx
           flange.f ~ 0 #TODO: model flow force
           connect(valve.port_a, port_a)
           connect(valve.port_b, port_b)
           valve.input.u ~ x * 2π * d]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [flange.v => 0])
end

@component function SpoolValve2Way(reversible = false; p_s_int, p_a_int, p_b_int, p_r_int,
                                   m, g, x_int, Cd, d, name)
    pars = @parameters begin
        p_s_int = p_s_int
        p_a_int = p_a_int
        p_b_int = p_b_int
        p_r_int = p_r_int

        m = m
        g = g

        x_int = x_int

        d = d

        Cd=Cd
    end

    vars = []

    systems = @named begin
        vSA = SpoolValve(reversible; p_a_int = p_s_int, p_b_int = p_a_int, x_int, Cd, d)
        vBR = SpoolValve(reversible; p_a_int = p_b_int, p_b_int = p_r_int, x_int, Cd, d)

        port_s = HydraulicPort(; p_int = p_s_int)
        port_a = HydraulicPort(; p_int = p_a_int)
        port_b = HydraulicPort(; p_int = p_b_int)
        port_r = HydraulicPort(; p_int = p_r_int)

        mass = Mass(; m = m, g = g)

        flange = MechanicalPort()
    end

    eqs = [connect(vSA.port_a, port_s)
           connect(vSA.port_b, port_a)
           connect(vBR.port_a, port_b)
           connect(vBR.port_b, port_r)
           connect(vSA.flange, vBR.flange, mass.flange, flange)]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [flange.v => 0])
end

@component function Actuator(; p_a_int, p_b_int, area_a, area_b, length_a_int, length_b_int,
                             m, g, x_int = 0, minimum_volume_a = 0, minimum_volume_b = 0,
                             name)
    pars = @parameters begin
        p_a_int = p_a_int
        p_b_int = p_b_int
        area_a = area_a
        area_b = area_b
        x_int = x_int
        length_a_int = length_a_int
        length_b_int = length_b_int
        minimum_volume_a = minimum_volume_a
        minimum_volume_b = minimum_volume_b
        m=m
        g=g
    end

    vars = @variables begin
        x(t) = x_int
        dx(t) = 0
    end

    systems = @named begin
        vol_a = DynamicVolume(; p_int = p_a_int, x_int = +x_int, area = area_a,
                              dead_volume = length_a_int * area_a,
                              minimum_volume = minimum_volume_a, direction = +1)
        vol_b = DynamicVolume(; p_int = p_b_int, x_int = -x_int, area = area_b,
                              dead_volume = length_b_int * area_b,
                              minimum_volume = minimum_volume_b, direction = -1)
        mass = Mass(; m, g)
        port_a = HydraulicPort(; p_int = p_a_int)
        port_b = HydraulicPort(; p_int = p_b_int)
        flange = MechanicalPort()
    end

    eqs = [connect(vol_a.port, port_a)
           connect(vol_b.port, port_b)
           connect(vol_a.flange, vol_b.flange, mass.flange, flange)
           x ~ vol_a.x
           dx ~ vol_a.dx]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [flange.v => 0])
end
