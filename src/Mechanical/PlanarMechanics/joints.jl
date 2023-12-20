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

https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Joints/Revolute.mo
"""
@component function Revolute(;
        name,
        constant_phi = nothing,
        constant_ω = nothing,
        constat_tau = nothing,
        use_flange = false)
    @named partial_frames = PartialTwoFrames()
    @unpack frame_a, frame_b = partial_frames
    @named fixed = Rotational.Fixed()
    systems = [frame_a, frame_b, fixed]

    vars = @variables begin
        phi(t) = 0.0
        ω(t) = 0.0
        α(t) = 0.0
        j(t) = 0.0
    end

    eqs = [
        phi ~ ifelse(constant_phi === nothing, phi, constant_phi),
        ω ~ ifelse(constant_ω === nothing, D(phi), constant_ω),
        α ~ ifelse(constant_ω === nothing, D(ω), 0.0),
        # rigidly connect positions
        frame_a.x ~ frame_b.x,
        frame_a.y ~ frame_b.y,
        frame_a.phi + phi ~ frame_b.phi,
        # balance forces
        frame_a.fx + frame_b.fx ~ 0,
        frame_a.fy + frame_b.fy ~ 0,
        # balance torques
        frame_a.j + frame_b.j ~ 0,
        j ~ ifelse(constat_tau === nothing, j, constat_tau),
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
  - `x`: [m] x-direction of the rod wrt. body system at phi=0
  - `y`: [m] y-direction of the rod wrt. body system at phi=0
  - `constant_f`: [N] Constant force in direction of elongation
  - `constant_s`: [m] Constant elongation of the joint"
  - `use_flange=false`: If `true`, a force flange is enabled, otherwise implicitly grounded"

# states
  - `s(t)`: [m] Elongation of the joint
  - `v(t)`: [m/s] Velocity of elongation
  - `a(t)`: [m/s²] Acceleration of elongation
  - `f(t)`: [N] Force in direction of elongation

# Connectors
  - `frame_a` [Frame](@ref)
  - `frame_b` [Frame](@ref)
  - `fixed` [Fixed](@ref) if `use_flange == false`
  - `flange_a` [Flange](@ref) if `use_flange == true`
  - `support` [Support](@ref) if `use_flange == true`

https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Joints/Prismatic.mo
"""
@component function Prismatic(;
        name,
        x,
        y,
        constant_f = 0,
        constant_s = 0,
        use_flange = false)
    @named partial_frames = PartialTwoFrames()
    @unpack frame_a, frame_b = partial_frames
    @named fixed = TranslationalModelica.Support()
    systems = [frame_a, frame_b, fixed]

    if use_flange
        @named flange_a = TranslationalModelica.Flange(; f = constant_f, constant_s)
        push!(systems, flange_a)
        @named support = TranslationalModelica.Support()
        push!(systems, support)
    end

    vars = @variables begin
        s(t) = 0.0
        v(t) = 0.0
        a(t) = 0.0
        f(t) = 0.0
    end

    R = [cos(frame_a.phi) -sin(frame_a.phi);
        sin(frame_a.phi) cos(frame_a.phi)]
    e0 = R * [x, y]
    r0 = e0 * s

    eqs = [
        # ifelse(constant_s === nothing, s ~ s, s ~ constant_s),
        ifelse(constant_f === nothing, f ~ f, f ~ constant_f),
        v ~ D(s),
        a ~ D(v),
        # rigidly connect positions
        frame_a.x + r0[1] ~ frame_b.x,
        frame_a.y + r0[2] ~ frame_b.y,
        frame_a.phi ~ frame_b.phi,
        frame_a.fx + frame_b.fx ~ 0,
        frame_a.fy + frame_b.fy ~ 0,
        frame_a.j + frame_b.j + r0[1] * frame_b.fy - r0[2] * frame_b.fx ~ 0,
        e0[1] * frame_a.fx + e0[2] * frame_a.fy ~ f,
    ]

    if use_flange
        push!(eqs, connect(fixed.flange, support))
    else
        # actutation torque
        push!(eqs, constant_f ~ 0)
    end

    pars = []

    return compose(ODESystem(eqs, t, vars, pars; name = name),
        systems...)
end
