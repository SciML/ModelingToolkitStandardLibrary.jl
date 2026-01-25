# Controlled Sources: VCVS, VCCS, CCVS, CCCS
# SPICE-equivalent controlled source elements for circuit analysis

"""
    VCVS(; name, G = 1.0)

Voltage Controlled Voltage Source (SPICE E-element).

A four-terminal device where the output voltage is proportional to the input voltage.
The input port has infinite impedance (no current flows into the sensing terminals).

# Parameters
- `G`: [V/V] Voltage gain (default: 1.0). Can be negative for inverting behavior.

# Connectors
- `p1` Positive input (sensing) terminal
- `n1` Negative input (sensing) terminal
- `p2` Positive output terminal
- `n2` Negative output terminal

# Equations
```
v2 = G × v1
i1 = 0  (infinite input impedance)
```

# Example
```julia
@named vcvs = VCVS(G = 10.0)  # Voltage amplifier with gain 10
```
"""
@component function VCVS(; name, G = 1.0)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    pars = @parameters begin
        G = G, [description = "Voltage gain [V/V]"]
    end

    eqs = Equation[
        v2 ~ G * v1,
        i1 ~ 0,
    ]

    sys = System(eqs, t, [], pars; name, systems = [])
    return extend(sys, twoport)
end

"""
    VCCS(; name, Gm = 0.001)

Voltage Controlled Current Source (SPICE G-element).

A four-terminal device where the output current is proportional to the input voltage.
The input port has infinite impedance (no current flows into the sensing terminals).

# Parameters
- `Gm`: [A/V or S] Transconductance (default: 0.001 S = 1 mS)

# Connectors
- `p1` Positive input (sensing) terminal
- `n1` Negative input (sensing) terminal
- `p2` Positive output terminal
- `n2` Negative output terminal

# Equations
```
i2 = Gm × v1
i1 = 0  (infinite input impedance)
```

# Example
```julia
@named vccs = VCCS(Gm = 0.01)  # Transconductance amplifier with Gm = 10 mS
```
"""
@component function VCCS(; name, Gm = 0.001)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    pars = @parameters begin
        Gm = Gm, [description = "Transconductance [A/V]"]
    end

    eqs = Equation[
        i2 ~ Gm * v1,
        i1 ~ 0,
    ]

    sys = System(eqs, t, [], pars; name, systems = [])
    return extend(sys, twoport)
end

"""
    CCVS(; name, Rm = 1000.0)

Current Controlled Voltage Source (SPICE H-element).

A four-terminal device where the output voltage is proportional to the input current.
The input port has zero impedance (no voltage drop across the sensing terminals).

# Parameters
- `Rm`: [V/A or Ω] Transimpedance (default: 1000.0 Ω = 1 kΩ)

# Connectors
- `p1` Positive input (sensing) terminal
- `n1` Negative input (sensing) terminal
- `p2` Positive output terminal
- `n2` Negative output terminal

# Equations
```
v2 = Rm × i1
v1 = 0  (zero input impedance)
```

# Example
```julia
@named ccvs = CCVS(Rm = 1000.0)  # Transimpedance amplifier with Rm = 1 kΩ
```
"""
@component function CCVS(; name, Rm = 1000.0)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    pars = @parameters begin
        Rm = Rm, [description = "Transimpedance [V/A]"]
    end

    eqs = Equation[
        v2 ~ Rm * i1,
        v1 ~ 0,
    ]

    sys = System(eqs, t, [], pars; name, systems = [])
    return extend(sys, twoport)
end

"""
    CCCS(; name, α = 1.0)

Current Controlled Current Source (SPICE F-element).

A four-terminal device where the output current is proportional to the input current.
The input port has zero impedance (no voltage drop across the sensing terminals).

# Parameters
- `α`: [A/A] Current gain (default: 1.0). Can be negative for inverting behavior.

# Connectors
- `p1` Positive input (sensing) terminal
- `n1` Negative input (sensing) terminal
- `p2` Positive output terminal
- `n2` Negative output terminal

# Equations
```
i2 = α × i1
v1 = 0  (zero input impedance)
```

# Example
```julia
@named cccs = CCCS(α = 100.0)  # Current amplifier with gain 100
```
"""
@component function CCCS(; name, α = 1.0)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    pars = @parameters begin
        α = α, [description = "Current gain [A/A]"]
    end

    eqs = Equation[
        i2 ~ α * i1,
        v1 ~ 0,
    ]

    sys = System(eqs, t, [], pars; name, systems = [])
    return extend(sys, twoport)
end
