"""
    PMSM(; name, R = 0.5, Ls = 1e-3, lambda = 0.01, p = 4, Jr = 1e-4, Dr = 1e-4)

Surface-mounted permanent-magnet synchronous machine (sinusoidal back-EMF,
isolated-wye stator, non-salient so `Ld = Lq = Ls`).

The three stator terminals are individual electrical [`Pin`](@ref)s; the star
point is internal and isolated (`i_a + i_b + i_c = 0`). The shaft is a
`Mechanical.Rotational.Flange` referenced to the housing
`Mechanical.Rotational.Support`, so the machine composes with the standard
`Mechanical.Rotational` inertia/load/sensor components like any other
electromechanical element (it is the magnet-excited, three-phase
generalisation of the `EMF` component).

Per-phase stator equation (`x ∈ {a, b, c}`, electrical angle
`θ = p·φ`, magnet flux linkage `ψ_x = lambda·cos(θ - kₓ·2π/3)`):

```
pin_x.v - v_star = R·i_x + Ls·d(i_x)/dt + d(ψ_x)/dt
```

Electromagnetic torque is taken in the flux form (no division by speed, so it
is well posed through standstill):

```
τ_e = -p·lambda·Σ i_x·sin(θ - kₓ·2π/3)
```

# States

  - `phi(t)`: [`rad`] Rotor mechanical angle relative to the support
  - `w(t)`: [`rad/s`] Rotor mechanical angular velocity
  - `i_a(t)`, `i_b(t)`, `i_c(t)`: [`A`] Stator phase currents
  - `tau_e(t)`: [`N·m`] Electromagnetic (air-gap) torque
  - `e_a(t)`, `e_b(t)`, `e_c(t)`: [`V`] Permanent-magnet back-EMF
  - `psi_a(t)`, `psi_b(t)`, `psi_c(t)`: [`Wb`] Stator flux linkage

# Connectors

  - `pin_a` Stator phase-A terminal ([`Pin`](@ref))
  - `pin_b` Stator phase-B terminal ([`Pin`](@ref))
  - `pin_c` Stator phase-C terminal ([`Pin`](@ref))
  - `flange` Rotor shaft (`Mechanical.Rotational.Flange`)
  - `support` Machine housing (`Mechanical.Rotational.Support`)

# Parameters

  - `R`: [`Ω`] Stator resistance per phase
  - `Ls`: [`H`] Synchronous (per-phase) inductance, `Ld = Lq = Ls`
  - `lambda`: [`Wb`] Permanent-magnet flux linkage (back-EMF constant)
  - `p`: Number of pole pairs
  - `Jr`: [`kg·m²`] Rotor moment of inertia
  - `Dr`: [`N·m·s`] Rotor viscous friction coefficient
"""
@component function PMSM(;
        name, R = 0.5, Ls = 1.0e-3, lambda = 0.01,
        p = 4, Jr = 1.0e-4, Dr = 1.0e-4
    )
    pars = @parameters begin
        R = R, [description = "Stator resistance per phase"]
        Ls = Ls, [description = "Synchronous (per-phase) inductance"]
        lambda = lambda, [description = "Permanent-magnet flux linkage"]
        p = p, [description = "Number of pole pairs"]
        Jr = Jr, [description = "Rotor moment of inertia"]
        Dr = Dr, [description = "Rotor viscous friction coefficient"]
    end

    systems = @named begin
        pin_a = Pin()
        pin_b = Pin()
        pin_c = Pin()
        flange = Flange()
        support = Support()
    end

    vars = @variables begin
        phi(t), [guess = 0.0, description = "Rotor mechanical angle"]
        w(t), [guess = 0.0, description = "Rotor mechanical speed"]
        i_a(t), [guess = 0.0, description = "Phase-A current"]
        i_b(t), [guess = 0.0, description = "Phase-B current"]
        i_c(t), [guess = 0.0, description = "Phase-C current"]
        v_star(t), [guess = 0.0, description = "Star-point potential"]
        tau_e(t), [guess = 0.0, description = "Electromagnetic torque"]
        e_a(t), [guess = 0.0, description = "Phase-A back-EMF"]
        e_b(t), [guess = 0.0, description = "Phase-B back-EMF"]
        e_c(t), [guess = 0.0, description = "Phase-C back-EMF"]
        psi_a(t), [guess = 0.0, description = "Phase-A stator flux linkage"]
        psi_b(t), [guess = 0.0, description = "Phase-B stator flux linkage"]
        psi_c(t), [guess = 0.0, description = "Phase-C stator flux linkage"]
    end

    g = 2 * pi / 3                       # 120° phase separation
    th = p * phi                         # electrical angle

    # Built as grouped sub-vectors so the explanatory comments stay on
    # plain statement lines (comments inside an array literal do not
    # survive autoformatting cleanly).

    # mechanical: shaft rigid to rotor, angle/speed referenced to support
    mech = [
        phi ~ flange.phi - support.phi,
        D(phi) ~ w,
    ]

    # stator phase currents flow in through the pins
    terminals = [
        i_a ~ pin_a.i,
        i_b ~ pin_b.i,
        i_c ~ pin_c.i,
    ]

    # per-phase stator voltage (phase-to-star); isolated wye closes v_star
    stator = [
        pin_a.v - v_star ~ R * i_a + Ls * D(i_a) + e_a,
        pin_b.v - v_star ~ R * i_b + Ls * D(i_b) + e_b,
        pin_c.v - v_star ~ R * i_c + Ls * D(i_c) + e_c,
        0 ~ i_a + i_b + i_c,
    ]

    # magnetic-domain observables: PM back-EMF e_x = d(ψ_pm,x)/dt and
    # total stator flux linkage ψ_x = Ls·i_x + lambda·cos(th - kₓ·g)
    magnetic = [
        e_a ~ -lambda * p * w * sin(th),
        e_b ~ -lambda * p * w * sin(th - g),
        e_c ~ -lambda * p * w * sin(th + g),
        psi_a ~ Ls * i_a + lambda * cos(th),
        psi_b ~ Ls * i_b + lambda * cos(th - g),
        psi_c ~ Ls * i_c + lambda * cos(th + g),
    ]

    # electromagnetic torque (flux form, no 1/w), rotor free body, and
    # housing reaction (massless stator) keeping Σ external τ = Jr·dw/dt
    dynamics = [
        tau_e ~ -p * lambda *
            (i_a * sin(th) + i_b * sin(th - g) + i_c * sin(th + g)),
        Jr * D(w) ~ tau_e - Dr * w + flange.tau,
        support.tau ~ tau_e - Dr * w,
    ]

    equations = [mech; terminals; stator; magnetic; dynamics]

    return System(equations, t, vars, pars; name, systems)
end
