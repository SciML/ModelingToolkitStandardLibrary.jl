@connector function HydraulicPort(; name, p_start = 1e5, m_flow_start = 0.0)
    @variables p(t) = p_start
    @variables m_flow(t) = m_flow_start [connect = Flow]
    return ODESystem(Equation[], t, [p, m_flow], []; name = name)
end
@doc """
    HydraulicPort(; name, p_start=1e5, m_flow_start=0.0)

Port for hydraulic systems.

# States:
- `p(t)`: [`Pa`] Pressure
- `m_flow(t)`: [`kg/s`] Mass flow rate

# Parameters:
- `p_start`: [`Pa`] Initial pressure
- `m_flow_start`: [`kg/s`] Initial mass flow rate
""" HydraulicPort

# smooth approximation of x * abs(x)
regAbs(x, delta = 0.01) = x * sqrt(x^2 + delta^2)

# Smoothed version of sign(x)*sqrt(abs(x)), which has a finite derivative at x=0
regRoot(x, delta = 0.01) = x / (x^2 + delta^2)^0.25

"""
    FluidProperties(;name, p_atm=101325.0, nu_atm=1.0034e-6, beta_atm=2.1791e9, rho_atm=998.21)

Defines the properties of the Liquid.

# Example:
```
@named fluid_props = FluidProperties(;rho_atm=998.21)
model = extend(model, fluid_props)
```

# Parameters:
- `p_atm = 101325.0`: [Pa] Atmospheric pressure
- `nu_atm = 1.0034e-6`: [m^2/s] Kinematic viscosity at atmospheric pressure
- `beta_atm = 2.1791e9`: [Pa] Isothermal bulk modulus at atmospheric pressure
- `rho_atm = 998.21`: [kg/m^3] Liquid density at atmospheric pressure
"""
function FluidProperties(;
    name,
    p_atm = 101325.0,
    nu_atm = 1.0034e-6,
    beta_atm = 2.1791e9,
    rho_atm = 998.21,
)
    pars = @parameters p_atm = p_atm nu_atm = nu_atm beta_atm = beta_atm rho_atm = rho_atm
    ODESystem(Equation[], t, [], pars; name = name)
end

"""
    calc_density(p, rho_atm, p_atm, beta_atm)

Computes the density of the liquid for the given conditions.

Gholizadeh, Hossein, Richard Burton, and Greg Schoenau. “Fluid Bulk Modulus: Comparison of Low 
Pressure Models.” International Journal of Fluid Power 13, no. 1 (January 2012): 7–16. 
https://doi.org/10.1080/14399776.2012.10781042.

Assumptions:
- Zero air in the liquid
- Constant bulk modulus

# Parameters:
- `rho_atm`: [kg/m^3] Liquid density at atmospheric pressure
- `p_atm`: [Pa] Atmospheric pressure
- `beta_atm`: [Pa] Isothermal bulk modulus at atmospheric pressure
"""
function calc_density(p, rho_atm, p_atm, beta_atm)
    return rho_atm * exp((p - p_atm) / beta_atm)
end
