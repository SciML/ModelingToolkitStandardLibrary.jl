"""
    NPN(;name, B_F, B_R, Is, V_T, V_A, phi_C, phi_E, Z_C, Z_E, Tau_f, Tau_r, C_jC0, C_jE0, C_CS, gamma_C, gamma_E, NF, NR)

Creates an NPN Bipolar Junction Transistor following a modified Ebers-Moll model. Includes an optional substrate pin and optional
Early voltage effect. 

    # Structural Parameters
        - `use_substrate`: If `true`, a substrate pin connector is available. If `false` it is 
        assumed the substrate is connected to the collector pin.

        - `use_Early`: If `true`, the Early effect is modeled, which takes in to account the effect 
        collector-base voltage variations have on the collector-base depletion region. In many cases this
        effectively means that the collector current has a dependency on the collector-emitter voltage.

        - `use_advanced_continuation`: When false, the `C_jC` and `C_jE` non-linear capacitance curves use 
        a simplified linear continuation starting when `V_BC` and `V_BE` are 0, respectively. If `true`, the `Z_C` and `Z_E` parameters 
        are used to start the linear continuation at `Phi_C - Z_C` and `Phi_E - Z_E`. 

    # Connectors
        - `b` Base Pin
        - `c` Collector Pin
        - `e` Emitter Pin
        - `s` Substrate Pin, only available when `use_substrate = true`

    # Parameters
        - `B_F`: Forward beta
        - `B_R`: Reverse beta
        - `Is`: Saturation current
        - `V_T`: Thermal voltage at 300K
        - `V_A`: Inverse Early voltage
        - `phi_C`: Collector junction exponent
        - `phi_E`: Emitter junction exponent
        - `Z_C`: Collector junction offset
        - `Z_E`: Emitter junction offset 
        - `Tau_f`: Forward transit time
        - `Tau_r`: Reverse transit time
        - `C_jC0`: Collector junction capacitance coefficient
        - `C_jE0`: Emitter junction capacitance coefficient
        - `C_CS`: Collector-substrate capacitance
        - `gamma_C`: Collector junction exponent
        - `gamma_E`: Emitter junction exponent
        - `NF`: Forward emission coefficient
        - `NR`: Reverse emission coefficient
"""
@component function NPN(; name,
                         use_substrate = false,
                         use_Early = true,
                         use_advanced_continuation = false,
                         V_BE = nothing, V_BC = nothing, ICC = nothing, IEC = nothing,
                         C_jC = nothing, C_jE = nothing, C_DC = nothing, C_DE = nothing,
                         I_sub = nothing, V_sub = nothing, V_CS = nothing,
                         B_F = 50.0, B_R = 0.1, Is = 1e-16, V_T = 0.026,
                         V_A = 0.02, phi_C = 0.8, phi_E = 0.6,
                         Z_C = 0.1, Z_E = 0.1,
                         Tau_f = 0.12e-9, Tau_r = 5e-9,
                         C_jC0 = 0.5e-12, C_jE0 = 0.4e-12,
                         C_CS = 1e-12,
                         gamma_C = 0.5, gamma_E = 1.0/3.0,
                         NF = 1.0, NR = 1.0)

    pars = @parameters begin
        B_F = B_F, [description = "Forward beta"]
        B_R = B_R, [description = "Reverse beta"]
        Is = Is, [description = "Saturation current"]
        V_T = V_T, [description = "Thermal voltage at 300K"]
        phi_C = phi_C, [description = "Collector junction scaling factor"]
        phi_E = phi_E, [description = "Emitter junction scaling factor"]
        Tau_f = Tau_f, [description = "Forward transit time"]
        Tau_r = Tau_r, [description = "Reverse transit time"]
        C_jC0 = C_jC0, [description = "Collector-junction capacitance coefficient"]
        C_jE0 = C_jE0, [description = "Emitter-junction capacitance coefficient"]
        C_CS = C_CS, [description = "Collector-substrate capacitance"]
        gamma_C = gamma_C, [description = "Collector junction exponent"]
        gamma_E = gamma_E, [description = "Emitter junction exponent"]
        NF = NF, [description = "Forward ideality exponent"]
        NR = NR, [description = "Reverse ideality exponent"]
    end

    if use_Early
        pars = [pars; @parameters(V_A = V_A, [description = "Inverse Early voltage"])]
    end

    if use_advanced_continuation
        pars = [pars; @parameters(Z_C = Z_C, [description = "Collector junction offset"])]
        pars = [pars; @parameters(Z_E = Z_E, [description = "Emitter junction offset"])]
    end

    @named b = Pin()
    @named e = Pin()
    @named c = Pin()
    systems = [b, e, c]

    if use_substrate
        @named s = Pin()
        push!(systems, s)
    end

    vars = @variables begin
        V_BE(t) = V_BE
        V_BC(t) = V_BC
        ICC(t) = ICC
        IEC(t) = IEC
        C_jC(t) = C_jC
        C_jE(t) = C_jE
        C_DC(t) = C_DC
        C_DE(t) = C_DE
        I_sub(t) = I_sub
        V_sub(t) = V_sub
        V_CS(t) = V_CS
    end

    equations = Equation[
        V_BE ~ b.v - e.v,
        V_BC ~ b.v - c.v,
        ICC ~ Is * (exp(V_BE / V_T) - 1),
        IEC ~ Is * (exp(V_BC / V_T) - 1)
    ]

    if !use_advanced_continuation
        push!(equations, C_jC ~ ifelse(V_BC / phi_C > 0.0, 1 + gamma_C * V_BC / phi_C,
            (C_jC0) / (1 - V_BC / phi_C)^gamma_C))
        push!(equations, C_jE ~ ifelse(V_BE / phi_E > 0.0, 1 + gamma_E * V_BE / phi_E,
            (C_jE0) / (1 - V_BE / phi_E)^gamma_E))
    end

    if use_advanced_continuation
        push!(equations, C_jC ~ ifelse(V_BC > phi_C - Z_C,
            ((C_jC0 * gamma_C * (1 - ((phi_C - Z_C) / phi_C))^(-gamma_C - 1)) / phi_C) *
            V_BC -
            ((C_jC0 * gamma_C * (1 - ((phi_C - Z_C) / phi_C))^(-gamma_C - 1)) / phi_C) *
            (phi_C - Z_C) + (C_jC0) / (1 - (phi_C - Z_C) / phi_C)^gamma_C,
            (C_jC0) / (1 - V_BC / phi_C)^gamma_C))

        push!(equations, C_jE ~ ifelse(V_BE > phi_E - Z_E,
            ((C_jE0 * gamma_E * (1 - ((phi_E - Z_E) / phi_E))^(-gamma_E - 1)) / phi_E) *
            V_BE -
            ((C_jE0 * gamma_E * (1 - ((phi_E - Z_E) / phi_E))^(-gamma_E - 1)) / phi_E) *
            (phi_E - Z_E) + (C_jE0) / (1 - (phi_E - Z_E) / phi_E)^gamma_E,
            (C_jE0) / (1 - V_BE / phi_E)^gamma_E))
    end

    push!(equations, C_DE ~ Tau_f * (Is / (NF * V_T)) * exp(V_BE / (NF * V_T)))
    push!(equations, C_DC ~ Tau_r * (Is / (NR * V_T)) * exp(V_BC / (NR * V_T)))

    if use_substrate
        push!(equations, s.i ~ I_sub)
        push!(equations, s.v ~ V_sub)
        push!(equations, V_CS ~ c.v - V_sub)
    end

    if !use_substrate
        push!(equations, V_sub ~ c.v)
    end

    push!(equations, I_sub ~ ifelse(use_substrate, -C_CS * D(V_CS), -C_CS * D(V_sub)))

    push!(equations, c.i ~ (ICC - IEC) * ifelse(use_Early, (1 - V_BC * V_A), 1.0) - IEC / B_R - (C_jC + C_DC) * D(V_BC) - I_sub)
    push!(equations, b.i ~ IEC / B_R + ICC / B_F + (C_jC + C_DC) * D(V_BC) + (C_jE + C_DE) * D(V_BE))
    push!(equations, e.i ~ -c.i - b.i - I_sub)

    return System(equations, t, vars, pars; name, systems)
end

"""
    PNP(;name, B_F, B_R, Is, V_T, V_A, phi_C, phi_E, Z_C, Z_E, Tau_f, Tau_r, C_jC0, C_jE0, C_CS, gamma_C, gamma_E, NF, NR)

Creates a PNP Bipolar Junction Transistor following a modified Ebers-Moll model. Includes an optional substrate pin and optional
Early voltage effect. 

    # Structural Parameters
        - `use_substrate`: If `true`, a substrate pin connector is available. If `false` it is 
        assumed the substrate is connected to the collector pin.

        - `use_Early`: If `true`, the Early effect is modeled, which takes in to account the effect 
        collector-base voltage variations have on the collector-base depletion region. In many cases this
        effectively means that the collector current has a dependency on the collector-emitter voltage.

        - `use_advanced_continuation`: When false, the `C_jC` and `C_jE` non-linear capacitance curves use 
        a simplified linear continuation starting when `V_CB` and `V_EB` are 0, respectively. If `true`, the `Z_C` and `Z_E` parameters 
        are used to start the linear continuation at `Phi_C - Z_C` and `Phi_E - Z_E`. 

    # Connectors
        - `b` Base Pin
        - `c` Collector Pin
        - `e` Emitter Pin
        - `s` Substrate Pin, only available when `use_substrate = true`

    # Parameters
        - `B_F`: Forward beta
        - `B_R`: Reverse beta
        - `Is`: Saturation current
        - `V_T`: Thermal voltage at 300K
        - `V_A`: Inverse Early voltage
        - `phi_C`: Collector junction exponent
        - `phi_E`: Emitter junction exponent
        - `Z_C`: Collector junction offset
        - `Z_E`: Emitter junction offset 
        - `Tau_f`: Forward transit time
        - `Tau_r`: Reverse transit time
        - `C_jC0`: Collector junction capacitance coefficient
        - `C_jE0`: Emitter junction capacitance coefficient
        - `C_CS`: Collector-substrate capacitance
        - `gamma_C`: Collector junction exponent
        - `gamma_E`: Emitter junction exponent
        - `NF`: Forward emission coefficient
        - `NR`: Reverse emission coefficient
"""
@component function PNP(; name,
                         use_substrate = false,
                         use_Early = true,
                         use_advanced_continuation = false,
                         V_EB = nothing, V_CB = nothing, ICC = nothing, IEC = nothing,
                         C_jC = nothing, C_jE = nothing, C_DC = nothing, C_DE = nothing,
                         I_sub = nothing, V_sub = nothing, V_CS = nothing,
                         B_F = 50.0, B_R = 0.1, Is = 1e-16, V_T = 0.026,
                         V_A = 0.02, phi_C = 0.8, phi_E = 0.6,
                         Z_C = 0.1, Z_E = 0.1,
                         Tau_f = 0.12e-9, Tau_r = 5e-9,
                         C_jC0 = 0.5e-12, C_jE0 = 0.4e-12,
                         C_CS = 1e-12,
                         gamma_C = 0.5, gamma_E = 1.0/3.0,
                         NF = 1.0, NR = 1.0)

    pars = @parameters begin
        B_F = B_F, [description = "Forward beta"]
        B_R = B_R, [description = "Reverse beta"]
        Is = Is, [description = "Saturation current"]
        V_T = V_T, [description = "Thermal voltage at 300K"]
        phi_C = phi_C, [description = "Collector junction scaling factor"]
        phi_E = phi_E, [description = "Emitter junction scaling factor"]
        Tau_f = Tau_f, [description = "Forward transit time"]
        Tau_r = Tau_r, [description = "Reverse transit time"]
        C_jC0 = C_jC0, [description = "Collector-junction capacitance coefficient"]
        C_jE0 = C_jE0, [description = "Emitter-junction capacitance coefficient"]
        C_CS = C_CS, [description = "Collector-substrate capacitance"]
        gamma_C = gamma_C, [description = "Collector junction exponent"]
        gamma_E = gamma_E, [description = "Emitter junction exponent"]
        NF = NF, [description = "Forward ideality exponent"]
        NR = NR, [description = "Reverse ideality exponent"]
    end

    if use_Early
        pars = [pars; @parameters(V_A = V_A, [description = "Inverse Early voltage"])]
    end

    if use_advanced_continuation
        pars = [pars; @parameters(Z_C = Z_C, [description = "Collector junction offset"])]
        pars = [pars; @parameters(Z_E = Z_E, [description = "Emitter junction offset"])]
    end

    @named b = Pin()
    @named e = Pin()
    @named c = Pin()
    systems = [b, e, c]

    if use_substrate
        @named s = Pin()
        push!(systems, s)
    end

    vars = @variables begin
        V_EB(t) = V_EB
        V_CB(t) = V_CB
        ICC(t) = ICC
        IEC(t) = IEC
        C_jC(t) = C_jC
        C_jE(t) = C_jE
        C_DC(t) = C_DC
        C_DE(t) = C_DE
        I_sub(t) = I_sub
        V_sub(t) = V_sub
        V_CS(t) = V_CS
    end

    equations = Equation[
        V_EB ~ e.v - b.v,
        V_CB ~ c.v - b.v,
        ICC ~ Is * (exp(V_EB / V_T) - 1),
        IEC ~ Is * (exp(V_CB / V_T) - 1)
    ]

    if !use_advanced_continuation
        push!(equations, C_jC ~ ifelse(V_CB / phi_C > 0.0, 1 + gamma_C * V_CB / phi_C,
            (C_jC0) / (1 - V_CB / phi_C)^gamma_C))
        push!(equations, C_jE ~ ifelse(V_EB / phi_E > 0.0, 1 + gamma_E * V_EB / phi_E,
            (C_jE0) / (1 - V_EB / phi_E)^gamma_E))
    end

    if use_advanced_continuation
        push!(equations, C_jC ~ ifelse(V_CB > phi_C - Z_C,
            ((C_jC0 * gamma_C * (1 - ((phi_C - Z_C) / phi_C))^(-gamma_C - 1)) / phi_C) *
            V_CB -
            ((C_jC0 * gamma_C * (1 - ((phi_C - Z_C) / phi_C))^(-gamma_C - 1)) / phi_C) *
            (phi_C - Z_C) + (C_jC0) / (1 - (phi_C - Z_C) / phi_C)^gamma_C,
            (C_jC0) / (1 - V_CB / phi_C)^gamma_C))

        push!(equations, C_jE ~ ifelse(V_EB > phi_E - Z_E,
            ((C_jE0 * gamma_E * (1 - ((phi_E - Z_E) / phi_E))^(-gamma_E - 1)) / phi_E) *
            V_EB -
            ((C_jE0 * gamma_E * (1 - ((phi_E - Z_E) / phi_E))^(-gamma_E - 1)) / phi_E) *
            (phi_E - Z_E) + (C_jE0) / (1 - (phi_E - Z_E) / phi_E)^gamma_E,
            (C_jE0) / (1 - V_EB / phi_E)^gamma_E))
    end

    push!(equations, C_DE ~ Tau_f * (Is / (NF * V_T)) * exp(V_EB / (NF * V_T)))
    push!(equations, C_DC ~ Tau_r * (Is / (NR * V_T)) * exp(V_CB / (NR * V_T)))

    if use_substrate
        push!(equations, s.i ~ I_sub)
        push!(equations, s.v ~ V_sub)
        push!(equations, V_CS ~ c.v - V_sub)
    end

    if !use_substrate
        push!(equations, V_sub ~ c.v)
    end

    push!(equations, I_sub ~ ifelse(use_substrate, -C_CS * D(V_CS), -C_CS * D(V_sub)))

    push!(equations, c.i ~ IEC / B_R - (ICC - IEC) * ifelse(use_Early, (1 - V_CB * V_A), 1.0) + (C_jC + C_DC) * D(V_CB) - I_sub)
    push!(equations, b.i ~ -IEC / B_R - ICC / B_F - (C_jC + C_DC) * D(V_CB) - (C_jE + C_DE) * D(V_EB))
    push!(equations, e.i ~ -c.i - b.i - I_sub)

    return System(equations, t, vars, pars; name, systems)
end
