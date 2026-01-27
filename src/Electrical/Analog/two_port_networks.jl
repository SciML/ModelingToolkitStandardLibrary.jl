# Two-Port Networks: IdealTransformer, Gyrator
# Power-conserving two-port network elements

"""
    IdealTransformer(; name, n = 1.0)

Ideal transformer with configurable turns ratio.

A lossless four-terminal device that transforms voltage and current between
its two ports according to the turns ratio. Power is conserved (P1 = -P2).

# Parameters
- `n`: [-] Turns ratio N1/N2 (default: 1.0). Can be any non-zero real number.

# Connectors
- `p1` Primary positive terminal
- `n1` Primary negative terminal
- `p2` Secondary positive terminal
- `n2` Secondary negative terminal

# Equations
```
v1 = n × v2         # Voltage transformation
n × i1 = -i2        # Current transformation (power conserving)
```

# Power Conservation
P1 = v1 × i1 = (n × v2) × (-i2/n) = -v2 × i2 = -P2

# Example
```julia
@named transformer = IdealTransformer(n = 10.0)  # 10:1 step-down transformer
```
"""
@component function IdealTransformer(; name, n = 1.0)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    pars = @parameters begin
        n = n, [description = "Turns ratio N1/N2"]
    end

    eqs = Equation[
        v1 ~ n * v2,
        n * i1 ~ -i2,
    ]

    sys = System(eqs, t, [], pars; name, systems = [])
    return extend(sys, twoport)
end

"""
    Gyrator(; name, R = 1000.0)

Ideal gyrator with configurable gyration resistance.

A lossless four-terminal device that performs impedance inversion. A gyrator
connected to a capacitor C presents an inductance L = R² × C at its input port.
This is useful for simulating inductors in integrated circuits.

# Parameters
- `R`: [Ω] Gyration resistance (default: 1000.0)

# Connectors
- `p1` Port 1 positive terminal
- `n1` Port 1 negative terminal
- `p2` Port 2 positive terminal
- `n2` Port 2 negative terminal

# Equations
```
v1 = R × i2
v2 = -R × i1
```

# Impedance Transformation
If port 2 is connected to impedance Z_load:
Z_in = R² / Z_load

For a capacitor C: Z_load = 1/(jωC)
Z_in = R² × jωC = jωL where L = R² × C

# Example
```julia
@named gyrator = Gyrator(R = 1000.0)  # 1 kΩ gyration resistance
# With 1μF capacitor: equivalent to 1H inductor
```
"""
@component function Gyrator(; name, R = 1000.0)
    @named twoport = TwoPort()
    @unpack v1, v2, i1, i2 = twoport

    pars = @parameters begin
        R = R, [description = "Gyration resistance [Ω]"]
    end

    eqs = Equation[
        v1 ~ R * i2,
        v2 ~ -R * i1,
    ]

    sys = System(eqs, t, [], pars; name, systems = [])
    return extend(sys, twoport)
end
