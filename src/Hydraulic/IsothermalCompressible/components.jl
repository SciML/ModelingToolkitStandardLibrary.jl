
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
    TubeBase(; p_int, area, length_int, effective_length_multiplier = 1, perimeter = 2 * sqrt(area * pi), shape_factor = 64, fluid_inertia_factor = 0, name)

Internal flow model of the fully developed flow friction, ignoring any compressibility.  Includes fluid inertia component to model wave propogation.  

# Parameters:
- `p_int`: [Pa] initial pressure
- `area`: [m^2] tube cross sectional area
- `length`: [m] real length of the tube
- `perimeter`: [m] perimeter of the pipe cross section (needed only for non-circular pipes)
- `Φ`: shape factor, see `friction_factor` function (set by optional `shape_factor` argument, needed only for non-circular pipes).  
- `Ε`: effective length multiplier, used to account for addition friction from flow development region and additional friction such as pipe bends, entrance/exit lossses, etc. (set by `effective_length_multiplier` argument)
- `fluid_inertia_factor`: account for wave propogation over the tube length, factor applied to mass flow derivative term

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
"""
@component function TubeBase(; p_int, area, length_int, effective_length_multiplier = 1, perimeter = 2 * sqrt(area * pi),
                             shape_factor = 64, fluid_inertia_factor = 0, name)
    pars = @parameters begin
        p_int = p_int
        area = area
        length_int=length_int
        perimeter = perimeter
        Φ = shape_factor
        Ε = effective_length_multiplier
        fluid_inertia_factor = fluid_inertia_factor
    end

    vars = @variables begin
        length(t) = length_int
        ddm(t) = 0
    end 

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

    eqs = [
        D(dm) ~ ddm
        Δp ~ 1 / 2 * ρ * u^2 * f * (length*Ε / d_h) + (length/area)*ddm*fluid_inertia_factor
        0 ~ port_a.dm + port_b.dm
        ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    Tube(N; p_int, area, length, effective_length=length, perimeter = 2 * sqrt(area * pi), shape_factor = 64, name)

Tube modeled with `N` segements which models the fully developed flow friction and compressibility (requires `N>1`).  When `N>1` the tube is segmented with `N` volumes and `N-1` resistive tube elements.

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
                         perimeter = 2 * sqrt(area * pi), fluid_inertia_factor = 0,
                         shape_factor = 64, name)
    @assert(N>0,
            "the Tube component must be defined with at least 1 segment (i.e. N>0), found N=$N")

    #TODO: How to set an assert effective_length >= length ??
    pars = @parameters begin
        p_int = p_int
        area = area
        length = length
        effective_length = effective_length
        perimeter = perimeter
        Φ = shape_factor
        fluid_inertia_factor = fluid_inertia_factor
    end

    vars = []

    ports = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
    end

    if N == 1

        @named pipe_base = TubeBase(; shape_factor=Φ, p_int, area, length_int=length, effective_length_multiplier=effective_length/length, fluid_inertia_factor, perimeter)

        eqs = [
            connect(pipe_base.port_a, port_a)
            connect(pipe_base.port_b, port_b)
            pipe_base.length ~ effective_length
            ]

        return ODESystem(eqs, t, vars, pars; name, systems = [ports; pipe_base])
    else
        pipe_bases = []
        for i in 1:(N - 1)
            x = TubeBase(; name = Symbol("p$i"), shape_factor = ParentScope(Φ),
                         p_int = ParentScope(p_int), area = ParentScope(area),
                         length_int = ParentScope(length) / (N - 1),
                         effective_length_multiplier = ParentScope(effective_length)/ParentScope(length),
                         fluid_inertia_factor = ParentScope(fluid_inertia_factor),
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
            push!(eqs, connect(volumes[i].port, pipe_bases[i - 1].port_b, pipe_bases[i].port_a)) 
        end
    
        for i = 1:(N-1)
            push!(eqs, pipe_bases[i].length ~ effective_length / (N - 1) )
        end
    
        return ODESystem(eqs, t, vars, pars; name, systems = [ports; pipe_bases; volumes])
    end

    
end
@deprecate Pipe Tube

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

    eqs = [connect(port_a, port_b, open.port)
           dm_a ~ port_a.dm
           dm_b ~ dm_a / n
           open.dm ~ dm_a - dm_b # extra flow dumps into an open port
           # port_b.dm ~ dm_b # divided flow goes to port_b
           ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

@component function ValveBase(reversible = false, directional = false; p_a_int, p_b_int,
                              area_int, Cd, name)
    pars = @parameters begin
        p_a_int = p_a_int
        p_b_int = p_b_int
        area_int = area_int
        Cd = Cd
    end

    systems = @named begin
        port_a = HydraulicPort(; p_int = p_a_int)
        port_b = HydraulicPort(; p_int = p_b_int)
    end

    vars = @variables area(t) = area_int

    # let
    ρ = full_density(port_a) # (full_density(port_a) + full_density(port_b)) / 2
    Δp = port_a.p - port_b.p
    dm = port_a.dm

    x = if reversible
        area
    else
        ifelse(area > 0, area, 0)
    end

    eqs = [
        0 ~ port_a.dm + port_b.dm,
    ]

    if directional
        push!(eqs, 0 ~ ifelse(Δp > 0, (dm) - (sqrt(2 * abs(Δp) * abs(ρ) / Cd) * x), Δp))
    else
        push!(eqs, sign(Δp) * dm ~ sqrt(2 * abs(Δp) * abs(ρ) / Cd) * x)
    end

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    Valve(reversible = false, directional=false; p_a_int, p_b_int, area_int, Cd, name)

Valve with `area` input and discharge coefficient `Cd` defined by https://en.wikipedia.org/wiki/Discharge_coefficient.  The input `directional` allows for 2 way flow restriction when `false`, and only applies restriction from `port_a` to `port_b` when true, making it like a check valve.

# Parameters:
- `p_a_int`: [Pa] initial pressure for `port_a`
- `p_b_int`: [Pa] initial pressure for `port_b`
- `area_int`: [m^2] initial valve opening
- `Cd`: discharge coefficient 

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
- `area`: real input setting the valve `area`.  When `reversible = true`, negative input reverses flow direction, otherwise a floor of 0 is enforced.
"""
@component function Valve(reversible = false, directional = false; p_a_int, p_b_int,
                          area_int, Cd, name)
    pars = @parameters begin
        p_a_int = p_a_int
        p_b_int = p_b_int
        area_int = area_int
        Cd = Cd
    end

    systems = @named begin
        port_a = HydraulicPort(; p_int = p_a_int)
        port_b = HydraulicPort(; p_int = p_b_int)
        area = RealInput()
        base = ValveBase(reversible; p_a_int, p_b_int, area_int, Cd)
    end

    vars = []

    eqs = [connect(base.port_a, port_a)
           connect(base.port_b, port_b)
           base.area ~ area.u]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [area.u => area_int])
end

@component function VolumeBase(; p_int, x_int = 0, area, dead_volume = 0, name)
    pars = @parameters begin
        p_int = p_int
        x_int = x_int
        area = area
        dead_volume = dead_volume
    end

    systems = @named begin port = HydraulicPort(; p_int) end

    vars = @variables begin
        x(t) = x_int
        dx(t) = 0
        rho(t) = liquid_density(port)
        drho(t) = 0
        vol(t) = dead_volume + area * x_int
    end

    # let
    dm = port.dm
    p = port.p

    eqs = [
           vol ~ dead_volume + area * x
           D(x) ~ dx
           D(rho) ~ drho
           rho ~ full_density(port, p)
           dm ~ drho * vol + rho * area * dx]

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

    systems = @named begin
        port = HydraulicPort(; p_int)
        vol = VolumeBase(; p_int, x_int = 0, area = 0.0, dead_volume = vol) #area set to 0.0, not needed or used here
    end

    eqs = [connect(vol.port, port)
           vol.dx ~ 0]

    ODESystem(eqs, t, [], pars; name, systems)
end




"""
    DynamicVolume(N, direction=+1; p_int, area, length, effective_length = length, fluid_inertia_factor = 0, perimeter = 2 * sqrt(area * pi), shape_factor = 64, minimum_volume = 0, damping_volume = 5 * minimum_volume, Cd = 1e4, name)

Volume with moving wall with `flange` connector for converting hydraulic energy to 1D mechanical.  The `direction` argument aligns the mechanical port with the hydraulic port, useful when connecting two dynamic volumes together in oppsing directions to create an actuator.  The `N` argument specifies the number of segments, with `N=1` the volume has equal pressure distribution, with `N>1` the volume is segmented by `N` and connected with `N-1` resistive tube elements.  The `DynamicVolume` also has a minimum volume feature with damping specified with the `minimum_volume` and `damping_volume` parameters.  When the minimum volume is reached the mass flow exiting the volume shuts off.  Mass flow can re-enter the volume without restriction.  A damping volume can be specified to provide a smooth linear transition from full flow to shut off.  

```
     ┌─────────────────┐ ───
     │                 │  ▲
                       │  │
dm ────►               │  │ area
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
- `vol(t)`: [m^3] volume `= x(t)*area + dead_volume`
- `p(t)`: [Pa] dynamic pressure

# Parameters:
- `p_int`: [Pa] initial pressure
- `area`: [m^2] moving wall area
- `perimeter`: [m] cross sectional perimeter, used to calculate the hydrualic diamter
- `length`: [m] starting length
- `effective_length`: [m] tube length to account for flow development and other restrictions
- `Φ`: shape factor, see `friction_factor` function (set by optional `shape_factor` argument, needed only for non-circular cylinders).  
- `fluid_inertia_factor`: account for wave propogation over the tube length, factor applied to mass flow derivative term
- `minimum_volume`: [m^3] if `vol(t) < minimum_volume` then mass flow `port.dm(t)` shuts off for exiting flow.
- `damping_volume`: [m^3] if `vol(t)` is between `damping_volume + minimum_volume` and `minimum_volume`, then a valve with restriction `Cd` closes linearly with area changing from 1 to 0 [m^2].  Restriction is applied to exiting flow only. 
- `Cd`: discharge coefficient (see Valve) for damping region


# Connectors:
- `port`: hydraulic port
- `flange`: mechanical translational port
"""
@component function DynamicVolume(N, direction=+1; p_int, area, length, effective_length = length, fluid_inertia_factor = 0,
    perimeter = 2 * sqrt(area * pi),
    shape_factor = 64, minimum_volume = 0, damping_volume = 5 * minimum_volume,
    Cd = 1e4, name)

    @assert(N>0, "the Tube component must be defined with more than 1 segment (i.e. N>1), found N=$N")
    @assert (direction == +1)||(direction == -1) "direction arument must be +/-1, found $direction"

    #TODO: How to set an assert effective_length >= length ??
    pars = @parameters begin
        p_int = p_int
        area = area
        length = length
        effective_length = effective_length
        perimeter = perimeter
        Φ = shape_factor
        fluid_inertia_factor = fluid_inertia_factor
        minimum_volume = minimum_volume
        damping_volume = damping_volume
        Cd = Cd
    end

    vars = @variables x(t)=0 vol(t)=length*area

    ports = @named begin
        port = HydraulicPort(; p_int)
        flange = MechanicalPort()
        damper = ValveBase(false, true; p_a_int = p_int, p_b_int = p_int, area_int = 1, Cd)
    end

    
    vol_ratio = (vol - minimum_volume) / damping_volume
    damper_area = ifelse(vol >= damping_volume + minimum_volume, 1,
                    ifelse((vol < damping_volume + minimum_volume) &
                        (vol > minimum_volume), vol_ratio, 0))  

    pipe_bases = []
    for i = 1:(N - 1)
        comp = TubeBase(; name = Symbol("p$i"), shape_factor = ParentScope(Φ),
                p_int = ParentScope(p_int), area = ParentScope(area),
                length_int = ParentScope(length) / (N - 1),
                effective_length_multiplier = ParentScope(effective_length)/ParentScope(length),
                fluid_inertia_factor = ParentScope(fluid_inertia_factor),
                perimeter = ParentScope(perimeter))
        push!(pipe_bases, comp)
    end

    #TODO: How to handle x_int?
    #TODO: Handle direction

    volumes = []
    for i = 1:N
        comp = VolumeBase(; name = Symbol("v$i"), p_int=ParentScope(p_int), x_int = 0, area=ParentScope(area), dead_volume = ParentScope(area)*ParentScope(length)/N)
        push!(volumes, comp)
    end

    eqs = [
        vol ~ (x + length)*area
        D(x) ~ flange.v*direction
        damper.area ~ damper_area
        connect(port, damper.port_b)
    ]

    if N == 1
        push!(eqs, connect(volumes[1].port, damper.port_a))
    else
        push!(eqs, connect(volumes[1].port, pipe_bases[1].port_a))
        push!(eqs, connect(volumes[end].port, pipe_bases[end].port_b, damper.port_a))
    end

    for i in 2:(N - 1)
        push!(eqs, connect(volumes[i].port, pipe_bases[i - 1].port_b, pipe_bases[i].port_a))
    end

    parts = []
    for i = 1:N
        push!(eqs, volumes[i].dx ~ ifelse( ( abs(x) >= (i-1)*(length/N) ) & ( abs(x) < (i)*(length/N) ) , +flange.v*direction, 0 ))

        push!(parts, ifelse( ( abs(x) >= (i-1)*(length/N) ) & ( abs(x) < (i)*(length/N) ) , (volumes[i].port.p)*(volumes[i].area),  0 ) )
        
    end

    for i=1:(N-1)
        push!(eqs, pipe_bases[i].length ~ volumes[i].vol/volumes[i].area)
    end

    push!(eqs, -flange.f*direction ~ +(parts...) ) 

    ODESystem(eqs, t, vars, pars; name, systems = [ports; pipe_bases; volumes])
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
        valve = ValveBase(reversible; p_a_int, p_b_int,
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
           valve.area ~ x * 2π * d]

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

        Cd = Cd
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

@component function Actuator(N; p_a_int, p_b_int, area_a, area_b, perimeter_a = 2 * sqrt(area_a * pi), perimeter_b = 2 * sqrt(area_b * pi), length_a_int, length_b_int,
                             shape_factor_a = 64, shape_factor_b = 64,
                             m, g, x_int = 0, minimum_volume_a = 0, minimum_volume_b = 0,
                             damping_volume_a = 0, damping_volume_b = 0,
                             fluid_inertia_factor = 0, 
                             Cd = 1e4,
                             name)
                             
    pars = @parameters begin
        p_a_int = p_a_int
        p_b_int = p_b_int
        area_a = area_a
        area_b = area_b
        perimeter_a = perimeter_a
        perimeter_b = perimeter_b
        shape_factor_a = shape_factor_a
        shape_factor_b = shape_factor_b
        x_int = x_int
        length_a_int = length_a_int
        length_b_int = length_b_int
        minimum_volume_a = minimum_volume_a
        minimum_volume_b = minimum_volume_b
        damping_volume_a = damping_volume_a
        damping_volume_b = damping_volume_b
        fluid_inertia_factor = fluid_inertia_factor
        Cd=Cd
        m = m
        g = g
    end

    vars = @variables begin
        x(t) = x_int
        dx(t) = 0
    end


    #TODO: include effective_length
    systems = @named begin
        vol_a = DynamicVolume(N, +1; p_int = p_a_int, area = area_a, perimeter = perimeter_a, shape_factor=shape_factor_a,
                              length = length_a_int,
                              minimum_volume = minimum_volume_a,
                              damping_volume = damping_volume_a, fluid_inertia_factor, Cd)
        vol_b = DynamicVolume(N, -1; p_int = p_b_int, area = area_b, perimeter = perimeter_b, shape_factor = shape_factor_b,
                              length = length_b_int,
                              minimum_volume = minimum_volume_b, perimeter = perimeter_b
                              damping_volume = damping_volume_b, fluid_inertia_factor, Cd)
        mass = Mass(; m, g)
        port_a = HydraulicPort(; p_int = p_a_int)
        port_b = HydraulicPort(; p_int = p_b_int)
        flange = MechanicalPort()
    end

    eqs = [connect(vol_a.port, port_a)
           connect(vol_b.port, port_b)
           connect(vol_a.flange, vol_b.flange, mass.flange, flange)
           x ~ vol_a.x
           dx ~ vol_a.flange.v]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [flange.v => 0])
end
