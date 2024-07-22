
"""
    Cap(; p_int, name)

Caps a hydraulic port to prevent mass flow in or out.

# Parameters:
- `p_int`: [Pa] initial pressure (set by `p_int` argument)

# Connectors:
- `port`: hydraulic port
"""
@component function Cap(; p_int, name)
    pars = @parameters p_int = p_int

    vars = @variables p(t) = p_int

    systems = @named begin
        port = HydraulicPort(; p_int = p_int)
    end

    eqs = [port.p ~ p
           port.dm ~ 0]

    ODESystem(eqs, t, vars, pars; name, systems)
end

@component function Open(; p_int, name)
    pars = @parameters p_int = p_int

    vars = @variables begin
        p(t) = p_int
        dm(t)
    end

    systems = @named begin
        port = HydraulicPort(; p_int = p_int)
    end

    eqs = [port.p ~ p
           port.dm ~ dm]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    TubeBase(add_inertia = true; p_int, area, length_int, head_factor = 1, perimeter = 2 * sqrt(area * pi), shape_factor = 64, name)

Variable length internal flow model of the fully developed incompressible flow friction.  Includes optional inertia term when `add_inertia = true` to model wave propagation.  Hydraulic ports have equal flow but variable pressure.  Density is averaged over the pressures, used to calculated average flow velocity and flow friction.

# States:
- `x`: [m] length of the pipe
- `ddm`: [kg/s^2] Rate of change of mass flow rate in control volume.

# Parameters:
- `p_int`: [Pa] initial pressure
- `area`: [m^2] tube cross sectional area
- `length_int`: [m] initial tube length
- `perimeter`: [m] perimeter of the pipe cross section (needed only for non-circular pipes)
- `shape_factor`: shape factor, see `friction_factor` function
- `head_factor`: effective length multiplier, used to account for addition friction from flow development and additional friction such as pipe bends, entrance/exit lossses, etc.

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
"""
@component function TubeBase(add_inertia = true, variable_length = true; p_int, area,
        length_int, head_factor = 1,
        perimeter = 2 * sqrt(area * pi),
        shape_factor = 64, name)
    pars = @parameters begin
        p_int = p_int
        area = area
        length_int = length_int
        perimeter = perimeter
        shape_factor = shape_factor
        head_factor = head_factor
    end

    @variables begin
        x(t) = length_int
        ddm(t)
    end

    vars = []
    if variable_length
        push!(vars, x)
        c = x
    else
        c = length_int
    end
    add_inertia && push!(vars, ddm)

    systems = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
    end

    # let ----------------------
    Δp = port_a.p - port_b.p
    dm = port_a.dm

    d_h = 4 * area / perimeter

    # Opting for a more numerically stable constant density (use head factor to compensate if needed)
    ρ = density_ref(port_a)  # (full_density(port_a) + full_density(port_b)) / 2
    μ = viscosity(port_a)

    f = friction_factor(dm, area, d_h, μ, shape_factor)
    u = dm / (ρ * area)

    shear = (1 / 2) * ρ * regPow(u, 2) * f * head_factor * (c / d_h)
    inertia = if add_inertia
        (c / area) * ddm
    else
        0
    end

    eqs = [0 ~ port_a.dm + port_b.dm
           domain_connect(port_a, port_b)]

    if variable_length
        push!(eqs, Δp ~ ifelse(c > 0, shear + inertia, zero(c)))
    else
        push!(eqs, Δp ~ shear + inertia)
    end

    if add_inertia
        push!(eqs, D(dm) ~ ddm)
    end

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    Tube(N, add_inertia=true; p_int, area, length, head_factor=1, perimeter = 2 * sqrt(area * pi), shape_factor = 64, name)

Constant length internal flow model discretized by `N` (`FixedVolume`: `N`, `TubeBase`:`N-1`) which models the fully developed flow friction, compressibility (when `N>1`), and inertia effects when `add_inertia = true`.  See `TubeBase` and `FixedVolume` for more information.

# Parameters:
- `p_int`: [Pa] initial pressure
- `area`: [m^2] tube cross sectional area
- `length`: [m] real length of the tube
- `perimeter`: [m] perimeter of the pipe cross section (needed only for non-circular pipes)
- `shape_factor`: shape factor, see `friction_factor` function
- `head_factor`: effective length multiplier, used to account for addition friction from flow development and additional friction such as pipe bends, entrance/exit lossses, etc.

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
"""
@component function Tube(N, add_inertia = true; p_int, area, length, head_factor = 1,
        perimeter = 2 * sqrt(area * pi),
        shape_factor = 64, name)
    @assert(N>0,
        "the Tube component must be defined with at least 1 segment (i.e. N>0), found N=$N")

    if N == 1
        return TubeBase(add_inertia,
            false;
            shape_factor,
            p_int,
            area,
            length_int = length,
            head_factor,
            perimeter,
            name)
    end

    #TODO: How to set an assert effective_length >= length ??
    pars = @parameters begin
        p_int = p_int
        area = area
        length = length
        head_factor = head_factor
        perimeter = perimeter
        shape_factor = shape_factor
    end

    vars = []

    ports = @named begin
        port_a = HydraulicPort(; p_int)
        port_b = HydraulicPort(; p_int)
    end

    pipe_bases = []
    for i in 1:(N - 1)
        x = TubeBase(add_inertia; name = Symbol("p$i"),
            shape_factor = ParentScope(shape_factor),
            p_int = ParentScope(p_int), area = ParentScope(area),
            length_int = ParentScope(length) / (N - 1),
            head_factor = ParentScope(head_factor),
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
        push!(eqs,
            connect(volumes[i].port, pipe_bases[i - 1].port_b, pipe_bases[i].port_a))
    end

    for i in 1:(N - 1)
        push!(eqs, pipe_bases[i].x ~ length / (N - 1))
    end

    return ODESystem(eqs, t, vars, pars; name, systems = [ports; pipe_bases; volumes])
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
        dm_a(t)
        dm_b(t)
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

@component function ValveBase(reversible = false; p_a_int, p_b_int, minimum_area = 0,
        area_int, Cd, Cd_reverse = Cd, name)
    pars = @parameters begin
        p_a_int = p_a_int
        p_b_int = p_b_int
        area_int = area_int
        Cd = Cd
        Cd_reverse = Cd_reverse
        minimum_area = minimum_area
    end

    systems = @named begin
        port_a = HydraulicPort(; p_int = p_a_int)
        port_b = HydraulicPort(; p_int = p_b_int)
    end

    vars = @variables begin
        area(t) = area_int
        y(t) = area_int
    end

    # let
    # Opting for a more numerically stable constant density (use head factor to compensate if needed)
    ρ = density_ref(port_a) #(full_density(port_a) + full_density(port_b)) / 2

    x = if reversible
        area
    else
        ifelse(area > minimum_area, area, minimum_area)
    end

    # let ------
    Δp = port_a.p - port_b.p
    dm = port_a.dm
    c = if reversible
        Cd
    else
        ifelse(Δp > 0, Cd, Cd_reverse)
    end

    eqs = [0 ~ port_a.dm + port_b.dm
           domain_connect(port_a, port_b)
           dm ~ regRoot(2 * Δp * ρ / c) * x
           y ~ x]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    Valve(reversible = false; p_a_int, p_b_int, area_int, Cd, Cd_reverse = Cd, minimum_area = 0, name)

Valve with `area` input and discharge coefficient `Cd` defined by https://en.wikipedia.org/wiki/Discharge_coefficient.  The `Cd_reverse` parameter allows for directional flow restriction, making it possible to define a check valve.

# Parameters:
- `p_a_int`: [Pa] initial pressure for `port_a`
- `p_b_int`: [Pa] initial pressure for `port_b`
- `area_int`: [m^2] initial valve opening
- `Cd`: discharge coefficient flowing from `a → b`
- `Cd_reverse`: discharge coefficient flowing from `b → a`
- `minimum_area`: when `reversible = false` applies a forced minimum area

# Connectors:
- `port_a`: hydraulic port
- `port_b`: hydraulic port
- `area`: real input setting the valve `area`.  When `reversible = true`, negative input reverses flow direction, otherwise a floor of `minimum_area` is enforced.
"""
@component function Valve(reversible = false; p_a_int, p_b_int,
        area_int, Cd, Cd_reverse = Cd,
        minimum_area = 0,
        name)
    pars = @parameters begin
        p_a_int = p_a_int
        p_b_int = p_b_int
        area_int = area_int
        Cd = Cd
        Cd_reverse = Cd_reverse
        minimum_area = minimum_area
    end

    systems = @named begin
        port_a = HydraulicPort(; p_int = p_a_int)
        port_b = HydraulicPort(; p_int = p_b_int)
        area = RealInput()
        base = ValveBase(reversible; p_a_int, p_b_int, area_int, Cd, Cd_reverse,
            minimum_area)
    end

    vars = []

    eqs = [connect(base.port_a, port_a)
           connect(base.port_b, port_b)
           base.area ~ area.u]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [area.u => area_int])
end

@component function VolumeBase(; p_int, x_int = 0, area, dead_volume = 0, Χ1 = 1, Χ2 = 1,
        name)
    pars = @parameters begin
        p_int = p_int
        x_int = x_int
        area = area
        dead_volume = dead_volume
    end

    systems = @named begin
        port = HydraulicPort(; p_int)
    end

    vars = @variables begin
        x(t) = x_int
        dx(t)
        rho(t) = liquid_density(port)
        drho(t)
        vol(t) = dead_volume + area * x_int
    end

    # let
    dm = port.dm
    p = port.p

    eqs = [vol ~ dead_volume + area * x
           D(x) ~ dx
           D(rho) ~ drho
           rho ~ full_density(port, p)
           dm ~ drho * vol * Χ1 + rho * area * dx * Χ2]

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
        p_int = p_int
        vol = vol
    end

    systems = @named begin
        port = HydraulicPort(; p_int)
    end

    vars = @variables begin
        rho(t) = liquid_density(port)
        drho(t)
    end

    # let
    dm = port.dm
    p = port.p

    eqs = [D(rho) ~ drho
           rho ~ full_density(port, p)
           dm ~ drho * vol]

    ODESystem(eqs, t, vars, pars; name, systems)
end

"""
    Volume(; x, dx=0, p, drho=0, dm=0, area, direction = +1, name)

Volume with moving wall with `flange` connector for converting hydraulic energy to 1D mechanical.  The `direction` argument aligns the mechanical port with the hydraulic port, useful when connecting two dynamic volumes together in oppsing directions to create an actuator.

```
     ┌─────────────────┐ ───
     │                 │  ▲
                       │  │
dm ────►               │  │ area
                       │  │
     │                 │  ▼
     └─────────────────┤ ───
                       │
                       └─► x (= ∫ flange.v * direction)
```

# Parameters:
## volume
- `p`: [Pa] initial pressure
- `area`: [m^2] moving wall area
- `x`: [m] initial wall position
- `dx=0`: [m/s] initial wall velocity
- `drho=0`: [kg/m^3/s] initial density derivative
- `dm=0`: [kg/s] initial flow

- `direction`: [+/-1] applies the direction conversion from the `flange` to `x`

# Connectors:
- `port`: hydraulic port
- `flange`: mechanical translational port

See also [`FixedVolume`](@ref), [`DynamicVolume`](@ref)
"""
@component function Volume(;
        #initial conditions
        x,
        dx = 0,
        p,
        drho = 0,
        dm = 0,

        #parameters
        area,
        direction = +1, name)
    pars = @parameters begin
        area = area
    end

    vars = @variables begin
        x(t) = x
        dx(t) = dx
        p(t) = p
        f(t) = p * area
        rho(t)
        drho(t) = drho
        dm(t) = dm
    end

    systems = @named begin
        port = HydraulicPort(; p_int = p)
        flange = MechanicalPort(; f, v = dx)
    end

    eqs = [
           # connectors
           port.p ~ p
           port.dm ~ dm
           flange.v * direction ~ dx
           flange.f * direction ~ -f

           # differentials
           D(x) ~ dx
           D(rho) ~ drho

           # physics
           rho ~ liquid_density(port, p)
           f ~ p * area
           dm ~ drho * x * area + rho * dx * area]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [rho => liquid_density(port)])
end

"""
DynamicVolume(N, add_inertia=true; p_int,  area, x_int = 0, x_max, x_min = 0, x_damp = x_min, direction = +1, perimeter = 2 * sqrt(area * pi), shape_factor = 64, head_factor = 1, Cd = 1e2, Cd_reverse = Cd, name)

Volume with moving wall with `flange` connector for converting hydraulic energy to 1D mechanical.  The `direction` argument aligns the mechanical port with the hydraulic port, useful when connecting two dynamic volumes together in oppsing directions to create an actuator.

```
     ┌─────────────────┐ ───
     │                 │  ▲
                       │  │
dm ────►               │  │ area
                       │  │
     │                 │  ▼
     └─────────────────┤ ───
                       │
                       └─► x (= ∫ flange.v * direction)
```

# Features:
- volume discretization with flow resistance and inertia: use `N` to control number of volume and resistance elements.  Set `N=0` to turn off volume discretization. See `TubeBase` for more information about flow resistance.
- minimum volume flow shutoff with damping and directional resistance.  Use `reversible=false` when problem defines volume position `x` and solves for `dm` to prevent numerical instability.

# Parameters:
## volume
- `p_int`: [Pa] initial pressure
- `area`: [m^2] moving wall area
- `x_int`: [m] initial wall position
- `x_max`: [m] max wall position, needed for volume discretization to apply the correct volume sizing as a function of `x`
- `x_min`: [m] wall position that shuts off flow and prevents negative volume.
- `x_damp`: [m] wall position that initiates a linear damping region before reaching full flow shut off.  Helps provide a smooth end stop.

- `direction`: [+/-1] applies the direction conversion from the `flange` to `x`

## flow resistance
- `perimeter`: [m] perimeter of the cross section (needed only for non-circular volumes)
- `shape_factor`: shape factor, see `friction_factor` function
- `head_factor`: effective length multiplier, used to account for addition friction from flow development and additional friction such as pipe bends, entrance/exit lossses, etc.

## flow shut off and damping
- `Cd`: discharge coefficient for flow out of the volume.  *Note: area is 1m² when valve is fully open.  Ensure this does not induce unwanted flow resistance.*
- `Cd_reverse`: discharge coefficient for flow into the volume. Use a lower value to allow easy wall release, in some cases the wall can "stick".


# Connectors:
- `port`: hydraulic port
- `flange`: mechanical translational port
"""
@component function DynamicVolume(N, add_inertia = true, reversible = false;
        p_int,
        area,
        x_int = 0,
        x_max,
        x_min = 0,
        x_damp = x_min,
        direction = +1,

        # Tube
        perimeter = 2 * sqrt(area * pi),
        shape_factor = 64,
        head_factor = 1,

        # Valve
        Cd = 1e2,
        Cd_reverse = Cd,
        minimum_area = 0,
        name)
    @assert(N>=0,
        "the Tube component must be defined with 0 or more segments (i.e. N>=0), found N=$N")
    @assert (direction == +1)||(direction == -1) "direction argument must be +/-1, found $direction"

    #TODO: How to set an assert effective_length >= length ??
    pars = @parameters begin
        p_int = p_int
        area = area

        x_int = x_int
        x_max = x_max
        x_min = x_min
        x_damp = x_damp

        # direction = direction

        perimeter = perimeter
        shape_factor = shape_factor
        head_factor = head_factor

        Cd = Cd
        Cd_reverse = Cd_reverse
        minimum_area = minimum_area
    end

    vars = @variables x(t)=x_int vol(t)=x_int * area

    ports = @named begin
        port = HydraulicPort(; p_int)
        flange = MechanicalPort(; f = -direction * p_int * area)
        damper = ValveBase(reversible;
            p_a_int = p_int,
            p_b_int = p_int,
            area_int = 1,
            Cd,
            Cd_reverse,
            minimum_area)
    end

    pipe_bases = []
    for i in 1:N
        comp = TubeBase(add_inertia; name = Symbol("p$i"),
            shape_factor = ParentScope(shape_factor),
            p_int = ParentScope(p_int), area = ParentScope(area),
            length_int = 0, #set in equations
            head_factor = ParentScope(head_factor),
            perimeter = ParentScope(perimeter))
        push!(pipe_bases, comp)
    end

    #TODO: How to handle x_int?
    #TODO: Handle direction
    @named moving_volume = VolumeBase(;
        p_int,
        x_int = 0,
        area,
        dead_volume = N == 0 ? area * x_int : 0,
        Χ1 = N == 0 ? 1 : 0,
        Χ2 = 1)

    ratio = (x - x_min) / (x_damp - x_min)

    damper_area = if reversible
        one(x)
    else
        ifelse(x >= x_damp, one(x), ifelse((x < x_damp) & (x > x_min), ratio, zero(x)))
    end

    eqs = [vol ~ x * area
           D(x) ~ flange.v * direction
           damper.area ~ damper_area
           connect(port, damper.port_b)]

    volumes = []
    if N > 0
        Δx = ParentScope(x_max) / N
        x₀ = ParentScope(x_int)

        for i in 1:N
            length = ifelse(x₀ > Δx * i,
                Δx,
                ifelse(x₀ - Δx * (i - 1) > 0,
                    x₀ - Δx * (i - 1),
                    zero(Δx)))

            comp = VolumeBase(; name = Symbol("v$i"), p_int = ParentScope(p_int),
                x_int = 0,
                area = ParentScope(area),
                dead_volume = ParentScope(area) * length, Χ1 = 1, Χ2 = 0)

            push!(volumes, comp)
        end

        push!(eqs, connect(moving_volume.port, volumes[end].port, pipe_bases[end].port_a))
        push!(eqs, connect(pipe_bases[1].port_b, damper.port_a))
        for i in 1:(N - 1)
            push!(eqs,
                connect(volumes[i].port, pipe_bases[i + 1].port_b, pipe_bases[i].port_a))
        end

        for i in 1:N
            push!(eqs,
                volumes[i].dx ~ ifelse(
                    (vol >= (i - 1) * (x_max / N) * area) &
                    (vol < (i) * (x_max / N) * area),
                    flange.v * direction, 0))
            push!(eqs, pipe_bases[i].x ~ volumes[i].vol / volumes[i].area)
        end
    else
        push!(eqs, connect(moving_volume.port, damper.port_a))
    end

    push!(eqs, moving_volume.dx ~ flange.v * direction)
    push!(eqs, -moving_volume.port.p * area * direction ~ flange.f)

    ODESystem(eqs, t, vars, pars; name,
        systems = [ports; pipe_bases; volumes; moving_volume],
        defaults = [flange.v => 0])
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
        dx(t)
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

@component function Actuator(N, add_inertia = true, reversible = false;
        p_a_int,
        p_b_int,
        area_a,
        area_b,
        perimeter_a = 2 * sqrt(area_a * pi),
        perimeter_b = 2 * sqrt(area_b * pi),
        length_a_int,
        length_b_int,
        shape_factor_a = 64,
        shape_factor_b = 64,
        head_factor_a = 1,
        head_factor_b = 1,
        m,
        g,
        x_int = 0,
        minimum_volume_a = 0,
        minimum_volume_b = 0,
        damping_volume_a = minimum_volume_a,
        damping_volume_b = minimum_volume_b,
        Cd = 1e4,
        Cd_reverse = Cd,
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
        head_factor_a = head_factor_a
        head_factor_b = head_factor_b
        x_int = x_int
        length_a_int = length_a_int
        length_b_int = length_b_int
        minimum_volume_a = minimum_volume_a
        minimum_volume_b = minimum_volume_b
        damping_volume_a = damping_volume_a
        damping_volume_b = damping_volume_b
        Cd = Cd
        Cd_reverse = Cd_reverse
        m = m
        g = g
    end

    vars = @variables begin
        x(t) = x_int
        dx(t)
    end

    total_length = length_a_int + length_b_int

    #TODO: include effective_length
    systems = @named begin
        vol_a = DynamicVolume(N, add_inertia, reversible; direction = +1,
            p_int = p_a_int,
            area = area_a,
            x_int = length_a_int,
            x_max = total_length,
            x_min = minimum_volume_a / area_a,
            x_damp = damping_volume_a / area_a,
            perimeter = perimeter_a,
            shape_factor = shape_factor_a,
            head_factor = head_factor_a,
            Cd,
            Cd_reverse)

        vol_b = DynamicVolume(N, add_inertia, reversible; direction = -1,
            p_int = p_b_int,
            area = area_b,
            x_int = length_b_int,
            x_max = total_length,
            x_min = minimum_volume_b / area_b,
            x_damp = damping_volume_b / area_b,
            perimeter = perimeter_b,
            shape_factor = shape_factor_b,
            head_factor = head_factor_b,
            Cd,
            Cd_reverse)
        mass = Mass(; m, g)
        port_a = HydraulicPort(; p_int = p_a_int)
        port_b = HydraulicPort(; p_int = p_b_int)
        flange = MechanicalPort()
    end

    eqs = [connect(vol_a.port, port_a)
           connect(vol_b.port, port_b)
           connect(vol_a.flange, vol_b.flange, mass.flange, flange)
           D(x) ~ dx
           dx ~ vol_a.flange.v]

    ODESystem(eqs, t, vars, pars; name, systems, defaults = [flange.v => 0])
end
