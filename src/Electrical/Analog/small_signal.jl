# Small-Signal Semiconductor Models: BJT_SmallSignal, MOSFET_SmallSignal
# Linearized transistor models for AC analysis

"""
    BJT_SmallSignal(; name, gm = 0.04, r_pi = 2500.0, r_o = 100000.0,
                     C_pi = 0.0, C_mu = 0.0)

Small-signal BJT model for AC analysis.

A linearized three-terminal transistor model with transconductance, input resistance,
output resistance, and optional junction capacitances. Suitable for frequency-domain
analysis and small-signal amplifier design.

# Parameters
- `gm`: [S] Transconductance (default: 0.04)
- `r_pi`: [Ω] Input resistance, base-emitter (default: 2500.0)
- `r_o`: [Ω] Output resistance (default: 100000.0)
- `C_pi`: [F] Base-emitter capacitance (default: 0.0, optional)
- `C_mu`: [F] Base-collector capacitance, Miller (default: 0.0, optional)

# Connectors
- `b` Base terminal
- `c` Collector terminal
- `e` Emitter terminal

# Equations
```
v_be = b.v - e.v
v_ce = c.v - e.v
b.i = v_be / r_pi + C_pi × D(v_be) + C_mu × D(b.v - c.v)
c.i = gm × v_be + v_ce / r_o - C_mu × D(b.v - c.v)
e.i = -b.i - c.i   # KCL
```

# Voltage Gain (Common-Emitter)
A_v ≈ -gm × (r_o ∥ R_C) for resistive collector load R_C

# Example
```julia
@named bjt = BJT_SmallSignal(gm = 0.04, r_pi = 2500.0, r_o = 100000.0)
```
"""
@component function BJT_SmallSignal(; name, gm = 0.04, r_pi = 2500.0, r_o = 100000.0,
        C_pi = 0.0, C_mu = 0.0)
    @named b = Pin()
    @named c = Pin()
    @named e = Pin()

    pars = @parameters begin
        gm = gm, [description = "Transconductance [S]"]
        r_pi = r_pi, [description = "Input resistance [Ω]"]
        r_o = r_o, [description = "Output resistance [Ω]"]
        C_pi = C_pi, [description = "Base-emitter capacitance [F]"]
        C_mu = C_mu, [description = "Base-collector capacitance [F]"]
    end

    vars = @variables begin
        v_be(t) = 0.0
        v_bc(t) = 0.0
        v_ce(t) = 0.0
    end

    eqs = Equation[
        v_be ~ b.v - e.v,
        v_bc ~ b.v - c.v,
        v_ce ~ c.v - e.v,
        b.i ~ v_be / r_pi + C_pi * D(v_be) + C_mu * D(v_bc),
        c.i ~ gm * v_be + v_ce / r_o - C_mu * D(v_bc),
        e.i ~ -b.i - c.i,
    ]

    return System(eqs, t, vars, pars; name, systems = [b, c, e])
end

"""
    MOSFET_SmallSignal(; name, gm = 0.01, r_ds = 50000.0,
                        C_gs = 0.0, C_gd = 0.0)

Small-signal MOSFET model for AC analysis.

A linearized three-terminal transistor model with transconductance, drain-source
resistance, and optional gate capacitances. Suitable for frequency-domain analysis
and small-signal amplifier design.

# Parameters
- `gm`: [S] Transconductance (default: 0.01)
- `r_ds`: [Ω] Drain-source resistance (default: 50000.0)
- `C_gs`: [F] Gate-source capacitance (default: 0.0, optional)
- `C_gd`: [F] Gate-drain capacitance, Miller (default: 0.0, optional)

# Connectors
- `g` Gate terminal
- `d` Drain terminal
- `s` Source terminal

# Equations
```
v_gs = g.v - s.v
v_ds = d.v - s.v
g.i = C_gs × D(v_gs) + C_gd × D(g.v - d.v)   # Gate current (capacitive only)
d.i = gm × v_gs + v_ds / r_ds - C_gd × D(g.v - d.v)
s.i = -g.i - d.i   # KCL
```

# Voltage Gain (Common-Source)
A_v ≈ -gm × (r_ds ∥ R_D) for resistive drain load R_D

# Example
```julia
@named mosfet = MOSFET_SmallSignal(gm = 0.01, r_ds = 50000.0)
```
"""
@component function MOSFET_SmallSignal(; name, gm = 0.01, r_ds = 50000.0,
        C_gs = 0.0, C_gd = 0.0)
    @named g = Pin()
    @named d = Pin()
    @named s = Pin()

    pars = @parameters begin
        gm = gm, [description = "Transconductance [S]"]
        r_ds = r_ds, [description = "Drain-source resistance [Ω]"]
        C_gs = C_gs, [description = "Gate-source capacitance [F]"]
        C_gd = C_gd, [description = "Gate-drain capacitance [F]"]
    end

    vars = @variables begin
        v_gs(t) = 0.0
        v_gd(t) = 0.0
        v_ds(t) = 0.0
    end

    eqs = Equation[
        v_gs ~ g.v - s.v,
        v_gd ~ g.v - d.v,
        v_ds ~ d.v - s.v,
        g.i ~ C_gs * D(v_gs) + C_gd * D(v_gd),
        d.i ~ gm * v_gs + v_ds / r_ds - C_gd * D(v_gd),
        s.i ~ -g.i - d.i,
    ]

    return System(eqs, t, vars, pars; name, systems = [g, d, s])
end
