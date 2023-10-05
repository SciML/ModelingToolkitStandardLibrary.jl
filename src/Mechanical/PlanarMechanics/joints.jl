"""
    Revolute(; name, phi = 0.0, tau = 0.0, use_flange = false)
A revolute joint

# parameters
  - `use_flange=false`: If `true`, a force flange is enabled, otherwise implicitly grounded"
  - `phi`: [rad] Intial angular position for the flange
  - `tau`: [N.m] Inital Cut torque in the flange

# states
  - `phi(t)`: [rad] angular position
  - `ω(t)`: [rad/s] angular velocity
  - `α(t)`: [rad/s²] angular acceleration
  - `j(t)`: [N.m] torque

# Connectors
  - `frame_a` [Frame](@ref)
  - `frame_b` [Frame](@ref)
  - `fixed` [Fixed](@ref) if `use_flange == false`
  - `flange_a` [Flange](@ref) if `use_flange == true`
  - `support` [Support](@ref) if `use_flange == true`

"""
@component function Revolute(; name, phi = 0.0, ω = 0.0, tau = 0.0, use_flange = false)
    @named partial_frames = PartialTwoFrames()
    @unpack frame_a, frame_b = partial_frames
    @named fixed = Rotational.Fixed()
    systems = [frame_a, frame_b, fixed]

    vars = @variables begin
        phi(t) = phi
        ω(t) = ω
        α(t) = 0.0
        j(t) = tau
    end

    eqs = [
        ω ~ D(phi),
        α ~ D(ω),
        # rigidly connect positions
        frame_a.x ~ frame_b.x,
        frame_a.y ~ frame_b.y,
        frame_a.phi + phi ~ frame_b.phi,
        # balance forces
        frame_a.fx + frame_b.fx ~ 0,
        frame_a.fy + frame_b.fy ~ 0,
        # balance torques
        frame_a.j + frame_b.j ~ 0,
        frame_a.j ~ j,
    ]

    if use_flange
        @named flange_a = Rotational.Flange(; phi, tau)
        push!(systems, flange_a)
        @named support = Rotational.Support()
        push!(systems, support)
        push!(eqs, connect(fixed.flange, support))
    else
        # actutation torque
        push!(eqs, j ~ 0)
    end

    pars = []

    return compose(ODESystem(eqs, t, vars, pars; name = name),
        systems...)
end

"""
    Prismatic(; name, rx, ry, f, s = 0, use_flange = false)
A prismatic joint

# parameters
  - `rx`: [m] x-direction of the rod wrt. body system at phi=0
  - `ry`: [m] y-direction of the rod wrt. body system at phi=0
  - `ex`: [m] x-component of unit vector in direction of r
  - `ey`: [m] y-component of unit vector in direction of r
  - `f`: [N] Force in direction of elongation
  - `s`: [m] Elongation of the joint"
  - `use_flange=false`: If `true`, a force flange is enabled, otherwise implicitly grounded"

# states
  - `s(t)`: [m] Elongation of the joint
  - `v(t)`: [m/s] Velocity of elongation
  - `a(t)`: [m/s²] Acceleration of elongation
  - `e0x(t)`: [m] x-component of unit vector resolved w.r.t inertial frame
  - `e0y(t)`: [m] y-component of unit vector resolved w.r.t inertial frame
  - `r0x(t)`: [m] x-component of the rod resolved w.r.t to inertal frame
  - `r0y(t)`: [m] y-length of the rod resolved w.r.t to inertal frame
  - `cos_phi(t)`: [degree] cos(phi)
  - `sin_phi(t)`: [degree] sin(phi)
    

# Connectors
  - `frame_a` [Frame](@ref)
  - `frame_b` [Frame](@ref)
  - `fixed` [Fixed](@ref) if `use_flange == false`
  - `flange_a` [Flange](@ref) if `use_flange == true`
  - `support` [Support](@ref) if `use_flange == true`
"""
@component function Prismatic(; name, rx, ry, ex, ey, f = 0, s = 0, use_flange = false)
    @named partial_frames = PartialTwoFrames()
    @unpack frame_a, frame_b = partial_frames
    @named fixed = TranslationalModelica.Support()
    systems = [frame_a, frame_b, fixed]

    if use_flange
        @named flange_a = TranslationalModelica.Flange(; f, s)
        push!(systems, flange_a)
        @named support = TranslationalModelica.Support()
        push!(systems, support)
    end

    vars = @variables begin
        s(t) = 0.0
        v(t) = 0.0
        a(t) = 0.0
        e0x(t)
        e0y(t)
        r0x(t)
        r0y(t)
        cos_phi(t)
        sin_phi(t)
    end

    eqs = [
        v ~ D(s),
        a ~ D(v),
        # rigidly connect positions
        frame_a.x + rx ~ frame_b.x,
        frame_a.y + ry ~ frame_b.y,
        frame_a.phi ~ frame_b.phi,
        # balance forces
        frame_a.fx + frame_b.fx ~ 0,
        frame_a.fy + frame_b.fy ~ 0,
        # balance torques
        cos_phi ~ cos(frame_a.phi),
        sin_phi ~ sin(frame_a.phi),
        e0x ~ cos_phi * ex - sin_phi * ey,
        e0y ~ sin_phi * ex + cos_phi * ey,
        r0x ~ e0x * s,
        r0y ~ e0y * s,
        frame_a.j + frame_b.j + r0x * (frame_b.fy - frame_a.fy) - r0y * (frame_b.fx - frame_a.fx) ~ 0,
        frame_a.fx * e0y - frame_a.fy * e0x ~ f,
    ]

    if use_flange
        push!(eqs, connect(fixed.flange, support))
    else
        # actutation torque
        push!(eqs, f ~ 0)
    end

    pars = []

    return compose(ODESystem(eqs, t, vars, pars; name = name),
        systems...)
end
