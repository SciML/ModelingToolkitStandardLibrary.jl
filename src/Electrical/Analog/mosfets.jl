"""
    NMOS(;name, V_tn, R_DS, lambda)

Creates an N-type MOSFET transistor

    # Structural Parameters
        - `use_transconductance`: If `true` the parameter `k_n` needs to be provided, and is used in the calculation of the current
        through the transistor. Otherwise, `mu_n`, `C_ox`, `W`, and `L` need to be provided and are used to calculate the transconductance.

        - `use_channel_length_modulation`: If `true` the channel length modulation effect is taken in to account. In essence this gives
        the drain-source current has a small dependency on the drains-source voltage in the saturation region of operation.

    # Connectors
        - `d` Drain Pin
        - `g` Gate Pin
        - `s` Source Pin

    # Parameters
        - `mu_n`: Electron mobility
        - `C_ox`: Oxide capacitance (F/m^2)
        - `W`: Channel width (m)
        - `L`: Channel length
        - `k_n`: MOSFET transconductance parameter

Based on the MOSFET models in (Sedra, A. S., Smith, K. C., Carusone, T. C., & Gaudet, V. C. (2021). Microelectronic circuits (8th ed.). Oxford University Press.)
"""
@component function NMOS(;
        name,
        use_transconductance = true,
        V_GS = nothing, V_DS = nothing, V_OV = nothing,
        V_tn = 0.8, R_DS = 1.0e7, lambda = 0.04,
        mu_n = nothing, C_ox = nothing, W = nothing, L = nothing,
        k_n = 20.0e-3
    )
    pars = @parameters begin
        V_tn = V_tn, [description = "Threshold voltage (V)"]
        R_DS = R_DS, [description = "Drain to source resistance (Ω)"]
        lambda = lambda, [description = "Channel length modulation coefficient (V^(-1))"]
    end

    if !use_transconductance
        pars = [pars; @parameters(mu_n = mu_n, [description = "Electron mobility"])]
        pars = [pars; @parameters(C_ox = C_ox, [description = "Oxide capacitance (F/m^2)"])]
        pars = [pars; @parameters(W = W, [description = "Channel width (m)"])]
        pars = [pars; @parameters(L = L, [description = "Channel length (m)"])]
    else
        pars = [pars; @parameters(k_n = k_n, [description = "MOSFET transconductance parameter"])]
    end

    @named d = Pin()
    @named g = Pin()
    @named s = Pin()
    systems = [d, g, s]

    # Calculate k_n if using individual parameters
    if !use_transconductance
        k_n = mu_n * C_ox * (W / L)
    end

    vars = @variables begin
        V_GS(t) = V_GS
        V_DS(t) = V_DS
        V_OV(t) = V_OV
    end

    equations = Equation[
        V_DS ~ ifelse(d.v < s.v, s.v - d.v, d.v - s.v),
        V_GS ~ g.v - ifelse(d.v < s.v, d.v, s.v),
        V_OV ~ V_GS - V_tn,
        d.i ~
            ifelse(d.v < s.v, -1, 1) * ifelse(
            V_GS < V_tn,
            V_DS / R_DS,
            ifelse(
                V_DS < V_OV,
                k_n * (1 + lambda * V_DS) * (V_OV - V_DS / 2) * V_DS + V_DS / R_DS,
                ((k_n * V_OV^2) / 2) * (1 + lambda * V_DS) + V_DS / R_DS
            )
        ),
        g.i ~ 0,
        s.i ~ -d.i,
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    PMOS(;name, V_tp, R_DS, lambda)

Creates an N-type MOSFET transistor

    # Structural Parameters
        - `use_transconductance`: If `true` the parameter `k_p` needs to be provided, and is used in the calculation of the current
        through the transistor. Otherwise, `mu_n`, `C_ox`, `W`, and `L` need to be provided and are used to calculate the transconductance.

        - `use_channel_length_modulation`: If `true` the channel length modulation effect is taken in to account. In essence this gives
        the drain-source current has a small dependency on the drains-source voltage in the saturation region of operation.

    # Connectors
        - `d` Drain Pin
        - `g` Gate Pin
        - `s` Source Pin

    # Parameters
        - `mu_p`: Electron mobility
        - `C_ox`: Oxide capacitance (F/m^2)
        - `W`: Channel width (m)
        - `L`: Channel length
        - `k_p`: MOSFET transconductance parameter

Based on the MOSFET models in (Sedra, A. S., Smith, K. C., Carusone, T. C., & Gaudet, V. C. (2021). Microelectronic circuits (8th ed.). Oxford University Press.)
"""
@component function PMOS(;
        name,
        use_transconductance = true,
        V_GS = nothing, V_DS = nothing,
        V_tp = -1.5, R_DS = 1.0e7, lambda = 1 / 25,
        mu_p = nothing, C_ox = nothing, W = nothing, L = nothing,
        k_p = 20.0e-3
    )
    pars = @parameters begin
        V_tp = V_tp, [description = "Threshold voltage (V)"]
        R_DS = R_DS, [description = "Drain-source resistance (Ω)"]
        lambda = lambda, [description = "Channel length modulation coefficient (V^(-1))"]
    end

    if !use_transconductance
        pars = [pars; @parameters(mu_p = mu_p, [description = "Hole mobility"])]
        pars = [pars; @parameters(C_ox = C_ox, [description = "Oxide capacitance (F/m^2)"])]
        pars = [pars; @parameters(W = W, [description = "Channel width (m)"])]
        pars = [pars; @parameters(L = L, [description = "Channel length (m)"])]
    else
        pars = [pars; @parameters(k_p = k_p, [description = "MOSFET transconductance parameter"])]
    end

    @named d = Pin()
    @named g = Pin()
    @named s = Pin()
    systems = [d, g, s]

    # Calculate k_p if using individual parameters
    if !use_transconductance
        k_p = mu_p * C_ox * (W / L)
    end

    vars = @variables begin
        V_GS(t) = V_GS
        V_DS(t) = V_DS
    end

    equations = Equation[
        V_DS ~ ifelse(d.v > s.v, s.v - d.v, d.v - s.v),
        V_GS ~ g.v - ifelse(d.v > s.v, d.v, s.v),
        d.i ~
            -ifelse(d.v > s.v, -1.0, 1.0) * ifelse(
            V_GS > V_tp,
            V_DS / R_DS,
            ifelse(
                V_DS > (V_GS - V_tp),
                k_p * (1 + lambda * V_DS) * ((V_GS - V_tp) - V_DS / 2) * V_DS +
                    V_DS / R_DS,
                ((k_p * (V_GS - V_tp)^2) / 2) * (1 + lambda * V_DS) + V_DS / R_DS
            )
        ),
        g.i ~ 0,
        s.i ~ -d.i,
    ]

    return System(equations, t, vars, pars; name, systems)
end
