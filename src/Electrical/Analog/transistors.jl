"""
    NPN(;name, B_F, B_R, Is, V_T, V_A, Phi_C, Phi_E, Z_C, Z_E, Tau_f, Tau_r, C_jC0, C_jE0, C_CS, gamma_C, gamma_E, NF, NR)

Creates an NPN Bipolar Junction Transistor following a modified Ebers-Moll model. 
"""



@mtkmodel NPN begin
    @variables begin
        V_BE(t)
        V_BC(t)
        ICC(t)
        IEC(t)

        C_jC(t)
        C_jE(t)
        C_DC(t)
        C_DE(t)

        I_sub(t)
        V_sub(t)
        V_CS(t)
    end

    @structural_parameters begin
        use_substrate = false
        use_Early = true
        use_advanced_continuation = false
    end

    @components begin
        b = Pin()
        e = Pin()
        c = Pin()

        if use_substrate
            s = Pin()
        end
    end

    @parameters begin
        B_F = 50.0
        B_R = 0.1
        Is = 1e-16
        V_T = 0.02

        if use_Early
            V_A = 0.02
        end

        Phi_C = 0.8
        Phi_E = 0.6

        if use_advanced_continuation
            Z_C = 0.1
            Z_E = 0.1
        end

        Tau_f = 0.12e-9
        Tau_r = 5e-9

        C_jC0 = 0.5e-12
        C_jE0 = 0.4e-12

        C_CS = 1e-12

        gamma_C = 0.5
        gamma_E = 1.0 / 3.0

        NF = 1.0
        NR = 1.0
    end

    @equations begin
        V_BE ~ b.v - e.v
        V_BC ~ b.v - c.v

        ICC ~ Is * (exp(V_BE / V_T) - 1)
        IEC ~ Is * (exp(V_BC / V_T) - 1)

        if !use_advanced_continuation
            C_jC ~ ifelse(V_BC / Phi_C > 0.0, 1 + gamma_C * V_BC / Phi_C,
                (C_jC0) / (1 - V_BC / Phi_C)^gamma_C)
            C_jE ~ ifelse(V_BE / Phi_E > 0.0, 1 + gamma_E * V_BE / Phi_E,
                (C_jE0) / (1 - V_BE / Phi_E)^gamma_E)
        end

        if use_advanced_continuation
            C_jC ~ if V_BC > Phi_C - Z_C
                ((C_jC0 * gamma_C * (1 - ((Phi_C - Z_C) / Phi_C))^(-gamma_C - 1)) / Phi_C) *
                V_BC -
                ((C_jC0 * gamma_C * (1 - ((Phi_C - Z_C) / Phi_C))^(-gamma_C - 1)) / Phi_C) *
                (Phi_C - Z_C) + (C_jC0) / (1 - (Phi_C - Z_C) / Phi_C)^gamma_C
            else
                (C_jC0) / (1 - V_BC / Phi_C)^gamma_C
            end

            C_jE ~ if V_BE > Phi_E - Z_E
                ((C_jE0 * gamma_E * (1 - ((Phi_E - Z_E) / Phi_E))^(-gamma_E - 1)) / Phi_E) *
                V_BE -
                ((C_jE0 * gamma_E * (1 - ((Phi_E - Z_E) / Phi_E))^(-gamma_E - 1)) / Phi_E) *
                (Phi_E - Z_E) + (C_jE0) / (1 - (Phi_E - Z_E) / Phi_E)^gamma_E
            else
                (C_jE0) / (1 - V_BE / Phi_E)^gamma_E
            end
        end

        C_DE ~ Tau_f * (Is / (NF * V_T)) * exp(V_BE / (NF * V_T))
        C_DC ~ Tau_r * (Is / (NR * V_T)) * exp(V_BC / (NR * V_T))

        if use_substrate
            s.i ~ I_sub
            s.v ~ V_sub
            V_CS ~ c.v - V_sub
        end

        if !use_substrate
            V_sub ~ c.v
        end

        I_sub ~ ifelse(use_substrate, -C_CS * D(V_CS), -C_CS * D(V_sub))

        C.i ~ (ICC - IEC) * ifelse(use_Early, (1 - V_BC * V_A), 1.0) - IEC / B_R -
              (C_jC + C_DC) * D(V_BC) - I_sub
        B.i ~ IEC / B_R + ICC / B_F + (C_jC + C_DC) * D(V_BC) + (C_jE + C_DE) * D(V_BE)
        E.i ~ -C.i - B.i - I_sub
    end
end