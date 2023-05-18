regPow(x, a, delta = 0.01) = x * (x * x + delta * delta)^((a - 1) / 2);
regRoot(x, delta = 0.01) = regPow(x, 0.5, delta)

"""
    HydraulicPort(;p_int, name)

Connector port for hydraulic components.

# Arguments:

- `p_int`: [Pa] initial gauge pressure

# States:
- `p`: [Pa] gauge total pressure
- `dm`: [kg/s] mass flow
"""
@connector function HydraulicPort(; p_int, name)
    pars = @parameters begin
        ρ
        β
        μ
        n
        ρ_gas
        p_gas
    end

    vars = @variables begin
        p(t) = p_int
        dm(t), [connect = Flow]
    end

    ODESystem(Equation[], t, vars, pars; name, defaults = [dm => 0])
end

"""
    HydraulicFluid(; density=997, bulk_modulus=2.09e9, viscosity=0.0010016, name)

Fluid parameter setter for isothermal compressible fluid domain.  Defaults given for water at 20°C and 0Pa gage (1atm absolute) reference pressure. Density is modeled using the Tait equation of state.  For pressures below the reference pressure, density is linearly interpolated to the gas state, this helps prevent pressures from going below the reference pressure.  

# Parameters:

- `ρ`: [kg/m^3] fluid density at 0Pa reference gage pressure (set by `density` argument)
- `Β`: [Pa] fluid bulk modulus describing the compressibility (set by `bulk_modulus` argument)
- `μ`: [Pa*s] or [kg/m-s] fluid dynamic viscosity  (set by `viscosity` argument)
- `n`: density exponent
- `ρ_gas`: [kg/m^3] density of fluid in gas state at reference gage pressure `p_gas` (set by `gas_density` argument)
- `p_gas`: [Pa] reference pressure (set by `gas_pressure` argument)
"""
@connector function HydraulicFluid(; density = 997, bulk_modulus = 2.09e9,
                                   viscosity = 0.0010016, gas_density = 0.0073955,
                                   gas_pressure = -1000, n = 1, name)
    pars = @parameters begin
        ρ = density
        β = bulk_modulus
        μ = viscosity
        n = n
        ρ_gas = gas_density
        p_gas = gas_pressure
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

    Re = density * u * d_h / viscosity
    f_laminar = shape_factor / Re

    Re = maximum([Re, 1])
    f_turbulent = (shape_factor / 64) * (0.79 * log(Re) - 1.64)^(-2)

    f = transition(2000, 3000, f_laminar, f_turbulent, Re)

    return f
end
@register_symbolic friction_factor(dm, area, d_h, density, viscosity, shape_factor)
Symbolics.derivative(::typeof(friction_factor), args, ::Val{1}) = 0
Symbolics.derivative(::typeof(friction_factor), args, ::Val{4}) = 0

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

density_ref(port) = port.ρ
gas_density_ref(port) = port.ρ_gas
gas_pressure_ref(port) = port.p_gas
bulk_modulus(port) = port.β
viscosity(port) = port.μ
liquid_density(port, p) = density_ref(port) * (1 + p / bulk_modulus(port))
liquid_density(port) = liquid_density(port, port.p)
function gas_density(port, p)
    density_ref(port) -
    p * (density_ref(port) - gas_density_ref(port)) / gas_pressure_ref(port)
end
full_density(port, p) = liquid_density(port, p)  #ifelse( p > 0, liquid_density(port, p), gas_density(port, p) )
full_density(port) = full_density(port, port.p)
