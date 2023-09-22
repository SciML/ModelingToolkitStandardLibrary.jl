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
@component function Revolute(; nane, phi = 0.0, tau = 0.0, use_flange = false)
    @named partial_frames = PartialTwoFrames()
    @unpack frame_a, frame_b = partial_frames
    @named fixed = Rotational.Fixed()

    if use_flange
        @named flange_a = Rotational.Flange(phi, tau)
        @named support = Rotational.Support()
    end

    vars = @varialbes begin
        phi = phi, [description = "Anugliar position"]
        ω = 0.0, [description = "Angular velocity"]
        α = 0.0, [description = "Angular acceleration"]
        j = tau, [description = "Torque"]
    end

    eqs = [
        ω ~ D(phi),
        α ~ D(ω),
        # rigidly connect position
        frame_a.x ~ frame_b.x,
        frame_a.y ~ frame_b.y,
        frame_a.phi + phi ~ frame_b.phi,
        # balance forces
        frame_a.fx + frame_b.fx = 0,
        frame_a.fy + frame_b.fy = 0,
        # balance torques
        frame_a.j + frame_b.j = 0,
        frame_a.j = j,
    ]

    if use_flange
        # actutation torque
        push!(eqs, connect(fixed.flange, support))
    else
        push!(eqs, j ~ ϕ)
    end

    pars = []

    return compose(ODESystem(eqs, t, vars, pars; name = name),
        partial_frames, fixed)
end