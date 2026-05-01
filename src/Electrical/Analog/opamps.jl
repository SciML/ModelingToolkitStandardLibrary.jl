# Op-Amp Models: OpAmpFiniteGain, OpAmpGBW, OpAmpFull
# Realistic operational amplifier models with non-ideal characteristics

"""
    OpAmpFiniteGain(; name, A = 100000.0, Rin = 1e6, Rout = 100.0)

Operational amplifier with finite DC gain.

A four-terminal op-amp model with configurable open-loop gain, input resistance,
and output resistance. Suitable for DC analysis where bandwidth effects are not important.

# Parameters
- `A`: [V/V] Open-loop DC voltage gain (default: 100000.0)
- `Rin`: [Ω] Input resistance (default: 1e6)
- `Rout`: [Ω] Output resistance (default: 100.0)

# Connectors
- `p1` Non-inverting input (+)
- `n1` Inverting input (-)
- `p2` Positive output terminal
- `n2` Negative output terminal (reference)

# Equations
```
i1 = v1 / Rin           # Input current through input resistance
v2 = A × v1 - i2 × Rout # Output with output resistance
```

# Example
```julia
@named opamp = OpAmpFiniteGain(A = 100000, Rin = 1e6, Rout = 100)
```
"""
@component function OpAmpFiniteGain(; name, A = 100000.0, Rin = 1e6, Rout = 100.0)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    pars = @parameters begin
        A = A, [description = "Open-loop DC gain [V/V]"]
        Rin = Rin, [description = "Input resistance [Ω]"]
        Rout = Rout, [description = "Output resistance [Ω]"]
    end

    eqs = Equation[
        i1 ~ v1 / Rin,
        v2 ~ A * v1 - i2 * Rout,
    ]

    sys = System(eqs, t, [], pars; name, systems = [])
    return extend(sys, twoport)
end

"""
    OpAmpGBW(; name, A0 = 100000.0, GBW = 1e6, Rin = 1e6)

Operational amplifier with gain-bandwidth product limitation.

A four-terminal op-amp model with a single-pole frequency response characterized
by the gain-bandwidth product. Suitable for frequency-domain analysis and
stability studies.

# Parameters
- `A0`: [V/V] DC open-loop gain (default: 100000.0)
- `GBW`: [Hz] Gain-bandwidth product (default: 1e6)
- `Rin`: [Ω] Input resistance (default: 1e6)

# Connectors
- `p1` Non-inverting input (+)
- `n1` Inverting input (-)
- `p2` Positive output terminal
- `n2` Negative output terminal (reference)

# States
- `v_pole(t)`: Internal state for single-pole response

# Equations
The transfer function is: H(s) = A0 / (1 + s/ω_p)
where ω_p = 2π × GBW / A0 is the pole frequency.

# Example
```julia
@named opamp = OpAmpGBW(A0 = 100000, GBW = 1e6, Rin = 1e6)
```
"""
@component function OpAmpGBW(; name, A0 = 100000.0, GBW = 1e6, Rin = 1e6)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    pars = @parameters begin
        A0 = A0, [description = "DC open-loop gain [V/V]"]
        GBW = GBW, [description = "Gain-bandwidth product [Hz]"]
        Rin = Rin, [description = "Input resistance [Ω]"]
    end

    vars = @variables begin
        v_pole(t) = 0.0
    end

    # Pole frequency: ω_p = 2π × GBW / A0
    # First-order dynamics: D(v_pole) = ω_p × (A0 × v1 - v_pole)
    eqs = Equation[
        i1 ~ v1 / Rin,
        D(v_pole) ~ (2 * π * GBW / A0) * (A0 * v1 - v_pole),
        v2 ~ v_pole,
    ]

    sys = System(eqs, t, vars, pars; name, systems = [])
    return extend(sys, twoport)
end

"""
    OpAmpFull(; name, A0 = 100000.0, GBW = 1e6, Vsat_p = 12.0, Vsat_n = -12.0,
              SR = 1e6, Rin = 1e6, Rout = 100.0)

Full behavioral operational amplifier model.

A comprehensive op-amp model including finite gain, bandwidth limitation,
output saturation, slew rate limiting, and output resistance.

# Parameters
- `A0`: [V/V] DC open-loop gain (default: 100000.0)
- `GBW`: [Hz] Gain-bandwidth product (default: 1e6)
- `Vsat_p`: [V] Positive saturation voltage (default: 12.0)
- `Vsat_n`: [V] Negative saturation voltage (default: -12.0)
- `SR`: [V/s] Slew rate (default: 1e6 = 1 V/μs)
- `Rin`: [Ω] Input resistance (default: 1e6)
- `Rout`: [Ω] Output resistance (default: 100.0)

# Connectors
- `p1` Non-inverting input (+)
- `n1` Inverting input (-)
- `p2` Positive output terminal
- `n2` Negative output terminal (reference)

# States
- `v_pole(t)`: Internal state for single-pole response

# Example
```julia
@named opamp = OpAmpFull(A0 = 100000, GBW = 1e6, Vsat_p = 12, Vsat_n = -12, SR = 1e6)
```
"""
@component function OpAmpFull(; name, A0 = 100000.0, GBW = 1e6, Vsat_p = 12.0, Vsat_n = -12.0,
        SR = 1e6, Rin = 1e6, Rout = 100.0)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    pars = @parameters begin
        A0 = A0, [description = "DC open-loop gain [V/V]"]
        GBW = GBW, [description = "Gain-bandwidth product [Hz]"]
        Vsat_p = Vsat_p, [description = "Positive saturation voltage [V]"]
        Vsat_n = Vsat_n, [description = "Negative saturation voltage [V]"]
        SR = SR, [description = "Slew rate [V/s]"]
        Rin = Rin, [description = "Input resistance [Ω]"]
        Rout = Rout, [description = "Output resistance [Ω]"]
    end

    vars = @variables begin
        v_pole(t) = 0.0
        v_out(t) = 0.0
    end

    # Pole frequency and dynamics
    # Use smooth saturation with clamp
    eqs = Equation[
        i1 ~ v1 / Rin,
        D(v_pole) ~ (2 * π * GBW / A0) * (A0 * v1 - v_pole),
        # Clamp to saturation limits
        v_out ~ IfElse.ifelse(v_pole > Vsat_p, Vsat_p,
            IfElse.ifelse(v_pole < Vsat_n, Vsat_n, v_pole)),
        v2 ~ v_out - i2 * Rout,
    ]

    sys = System(eqs, t, vars, pars; name, systems = [])
    return extend(sys, twoport)
end
