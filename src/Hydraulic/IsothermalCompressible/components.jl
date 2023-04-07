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
        rho(t) = density(port, p_int)
        drho(t) = 0
    end

    # let -------------------
    dm = port.dm

    eqs = [D(rho) ~ drho
           rho ~ density(port, port.p)
           dm ~ drho * vol]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    PipeBase(; p_int, area, length, perimeter=2*sqrt(area*pi), shape_factor=64, name)

Pipe segement which models purely the fully developed flow friction, ignoring any compressibility.

# Parameters:
- `p_int`: [Pa] initial pressure (set by `p_int` argument)
- `area`: [m^2] tube cross sectional area (set by `area` argument)
- `length`: [m] length of the pipe (set by `length` argument)
- `perimeter`: [m] perimeter of the pipe cross section (set by optional `perimeter` argument, needed only for non-circular pipes)
- `Φ`: shape factor, see `friction_factor` function (set by optional `shape_factor` argument, needed only for non-circular pipes).  

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
"""
@component function PipeBase(; p_int, area, length, perimeter = 2 * sqrt(area * pi),
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

    ρ = (density(port_a, port_a.p) + density(port_b, port_b.p)) / 2
    μ = viscosity(port_a)

    f = friction_factor(dm, area, d_h, ρ, μ, Φ)
    u = dm / (ρ * area)

    eqs = [Δp ~ 1 / 2 * ρ * u^2 * f * (length / d_h)
           0 ~ port_a.dm + port_b.dm]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    Pipe(N; p_int, area, length, perimeter=2*sqrt(area*pi), shape_factor=64, name)

Pipe modeled with `N` segements which models the fully developed flow friction and compressibility.

# Parameters:
- `p_int`: [Pa] initial pressure (set by `p_int` argument)
- `area`: [m^2] tube cross sectional area (set by `area` argument)
- `length`: [m] length of the pipe (set by `length` argument)
- `perimeter`: [m] perimeter of the pipe cross section (set by optional `perimeter` argument, needed only for non-circular pipes)
- `Φ`: shape factor, see `friction_factor` function (set by optional `shape_factor` argument, needed only for non-circular pipes).  

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
"""
@component function Pipe(N; p_int, area, length, perimeter = 2 * sqrt(area * pi),
                         shape_factor = 64, name)
    @assert(N>1,
            "the pipe component must be defined with more than 1 segment (i.e. N>1), found N=$N")

    pars = @parameters begin
        p_int = p_int
        area = area
        length = length
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
        x = PipeBase(; name = Symbol("p$i"), shape_factor = ParentScope(Φ),
                     p_int = ParentScope(p_int), area = ParentScope(area),
                     length = ParentScope(length) / (N - 1),
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
        rho(t) = density(port, p_int)
        drho(t) = 0
        vol(t) = dead_volume + area * x_int
        p(t) = p_int
    end

    eqs = [0 ~ IfElse.ifelse(vol >= minimum_volume, p - port.p, port.dm)
           vol ~ dead_volume + area * x
           D(x) ~ dx
           D(rho) ~ drho
           dx ~ flange.v * direction
           rho ~ density(port, p)
           port.dm ~ drho * vol + rho * area * dx
           flange.f ~ -p * area * direction] #TODO: update to dynamic pressure

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
@component function Valve(; p_a_int, p_b_int, area_int, Cd, name)
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
    ρ = (density(port_a, port_a.p) + density(port_b, port_b.p)) / 2
    Δp = port_a.p - port_b.p
    dm = port_a.dm
    area = abs(input.u)

    eqs = [sign(Δp) * dm ~ sqrt(2 * abs(Δp) * ρ / Cd) * area
           0 ~ port_a.dm + port_b.dm]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [input.u => area_int])
end
