"""
    HydraulicPort(;p_int, name)

Connector port for hydraulic components.

# Arguments:

- `p_int`: [Pa] initial gauge pressure

# States:
- `p`: [Pa] gauge pressure
- `dm`: [kg/s] mass flow
"""
@connector function HydraulicPort(; p_int, name)
    pars = @parameters begin
        ρ
        β
        μ
    end

    vars = @variables begin
        p(t) = p_int
        dm(t), [connect = Flow]
    end

    ODESystem(Equation[], t, vars, pars; name, defaults = [dm => 0])
end

"""
    HydraulicFluid(; density=997, bulk_modulus=2.09e9, viscosity=0.0010016, name)

Fluid parameter setter for isothermal compressible fluid domain.  Defaults given for water at 20°C and 1atm reference pressure.

# Parameters:

- `ρ`: [kg/m^3] fluid density at reference pressure (set by `density` argument)
- `Β`: [Pa] fluid bulk modulus describing the compressibility (set by `bulk_modulus` argument)
- `μ`: [Pa*s] or [kg/m-s] fluid dynamic viscosity  (set by `viscosity` argument)
"""
@connector function HydraulicFluid(; density = 997, bulk_modulus = 2.09e9,
                                   viscosity = 0.0010016, name)
    pars = @parameters begin
        ρ = density
        β = bulk_modulus
        μ = viscosity
    end

    vars = @variables begin dm(t), [connect = Flow] end

    eqs = [
        dm ~ 0,
    ]

    ODESystem(eqs, t, vars, pars; name, defaults = [dm => 0])
end

"""
    friction_factor(dm, area, d_h, density, viscosity, shape_factor)

Calculates the friction factor ``f`` for fully developed flow in a tube such that ``Δp = f \\cdot \\rho \\frac{u^2}{2} \\frac{l}{d_h}`` where 

- ``Δp``: [Pa] is the pressure difference over the tube length ``l``
- ``\\rho``: [kg/m^3] is the average fluid density
- ``u``: [m/s] is the average fluid velocity
- ``l``: [m] is the tube length 

The friction factor is calculated for laminar and turbulent flow with a transition region between Reynolds number 2000 to 3000.  Turbulent flow equation is for smooth tubes, valid for the Reynolds number range up to 5e6.

# Arguments:

- `dm`: [kg/s] mass flow
- `area`: [m^2] tube cross sectional area
- `d_h`: [m] tube hydraulic diameter.  For circular tubes d_h is the tube diameter, otherwise it can be found from `4*area/perimeter`
- `density`: [kg/m^3] fluid density
- `viscosity`: [Pa*s] or [kg/m-s] fluid dynamic viscosity
- `shape_factor`: the constant defining the laminar fully developed constant f*Re related to the shape of the tube cross section

Reference: Introduction to Fluid Mechanics, Fox & McDonald, 5th Edition, equations 8.19 and 8.21
"""
function friction_factor(dm, area, d_h, density, viscosity, shape_factor)
    u = abs(dm) / (density * area)
    Re = maximum([density * u * d_h / viscosity, 1])

    f_laminar = shape_factor / Re
    f_turbulent = (shape_factor / 64) * (0.79 * log(Re) - 1.64)^(-2)

    f = transition(2000, 3000, f_laminar, f_turbulent, Re)

    return f
end
@register_symbolic friction_factor(dm, area, d_h, density, viscosity, shape_factor)

function transition(x1, x2, y1, y2, x)
    if x < x1
        return y1
    elseif x > x2
        return y2
    else
        u = (x - x1) / (x2 - x1)
        blend = 3 * u^2 - 2 * u^3
        return (1 - blend) * y1 + blend * y2
    end
end

density(port) = port.ρ
bulk_modulus(port) = port.β
viscosity(port) = port.μ
density(port, p) = density(port) * (1 + p / bulk_modulus(port))
