"""
    Reservoir(; name, p)

Infinite reservoir; set BC for hydraulic system

# Connectors:
- `a` port [`HydraulicPort`](@ref)

# Parameters:
- `p`: [`Pa`] Pressure of the reservoir
"""
function Reservoir(; name, p)
    @parameters p = p
    @named a = HydraulicPort()
    eqs = [a.p ~ p]
    return ODESystem(eqs, t, [], [p]; systems = [a], name = name)
end

"""
    LocalRestriction(; name, A, Cd, Re_crit)

# States:

# Connectors:
- `a` left port [`HydraulicPort`](@ref)
- `b` right port [`HydraulicPort`](@ref)

# Parameters:
- `A`: [`m^2`] Effective area of the restriction
- `Cd`: [] Discharge coefficient
- `Re_crit=150`: [] Critical Reynolds number
"""
function LocalRestriction(; name, A, Cd, Re_crit = 150)
    @named a = HydraulicPort()
    @named b = HydraulicPort()

    pars = @parameters begin
        A = A
        Cd = Cd
        p_atm# , [scope=:parent]
        nu_atm# , [scope=:parent]
        beta_atm# , [scope=:parent]
        rho_atm# , [scope=:parent]
    end
    # lift the fluid parameters into parent scope
    p_atm = ParentScope(p_atm)
    nu_atm = ParentScope(nu_atm)
    beta_atm = ParentScope(beta_atm)
    rho_atm = ParentScope(rho_atm)

    pars = [A, Cd, p_atm, nu_atm, beta_atm, rho_atm]

    rho = calc_density((b.p + a.p) / 2, rho_atm, p_atm, beta_atm)

    Δp = b.p - a.p
    p_cr = pi / 4 * rho / (2 * A) * (Re_crit * nu_atm / Cd)^2

    eqs = [
        # Momentum balance
        b.m_flow ~ Cd * A * sqrt(2 * rho) * regRoot(Δp, p_cr)
        # Mass balance
        a.m_flow + b.m_flow ~ 0
    ]
    return ODESystem(eqs, t, [], pars; name = name, systems = [a, b])
end

"""
    ConstantVolume(; name, V, p_start)

Constant volume chamber.

# States:
- `p`: [`Pa`] Pressure inside the volume

# Connectors:
- `a` port [`HydraulicPort`](@ref)

# Parameters:
- `V`: [`m^3`] Constant volume
- `p_start=0.0`: [`Pa`] Initial pressure inside the volume
"""
function ConstantVolume(; name, V, p_start = 0.0)
    @named a = HydraulicPort()
    @variables p(t) = p_start # pressure inside the volume
    pars = @parameters begin
        V = V
        p_atm# , [scope=:parent]
        nu_atm# , [scope=:parent]
        beta_atm# , [scope=:parent]
        rho_atm# , [scope=:parent]
    end
    # lift the fluid parameters into parent scope
    p_atm = ParentScope(p_atm)
    nu_atm = ParentScope(nu_atm)
    beta_atm = ParentScope(beta_atm)
    rho_atm = ParentScope(rho_atm)

    pars = [V, p_atm, nu_atm, beta_atm, rho_atm]

    rho = calc_density(p, rho_atm, p_atm, beta_atm)

    eqs = [
        # Mass conservation
        D(p) ~ 1 / (Symbolics.derivative(rho, p) * V) * a.m_flow
        a.p ~ p # no pressure loss from inlet to volume
    ]

    return ODESystem(eqs, t, [p], pars; name = name, systems = [a])
end
