"""
    Fixed(; name, r = (0.0, 0.0), phi = 0.0)

Frame fixed in the planar world frame at a given position and orientation

# Parameters:

  - `x`: [m] Fixed absolute x-position, resolved in planarWorld frame
  - `y`: [m] Fixed absolute y-position, resolved in planarWorld frame
  - `phi`: [rad] Fixed angle

# Connectors:

  - `frame: 2-dim. Coordinate system

https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Parts/Fixed.mo
"""
@mtkmodel Fixed begin
    @parameters begin
        x = 0, [description = "Fixed absolute x-position, resolved in planarWorld frame"]
        y = 0, [description = "Fixed absolute y-position, resolved in planarWorld frame"]
        phi = 0, [description = "Fixed angle"]
    end

    @components begin
        frame = Frame()
    end

    @equations begin
        frame.x ~ x
        frame.y ~ y
        frame.phi ~ phi
    end
end

"""
    Body(; name, m, j, r, g = nothing)

Body component with mass and inertia

# Parameters:

  - `m`: [kg] mass of the body
  - `j`: [kg.m²] inertia of the body with respect to the origin of `frame` along the z-axis of `frame`
  - `r`: [m, m] (optional) Translational position x,y-position
  - `gy`: [m/s²] (optional) gravity field acting on the mass in the y-direction, positive value acts in the positive direction defaults to -9.807

# States:

  - `rx`: [m] x position
  - `ry`: [m] y position
  - `vx`: [m/s] x velocity
  - `vy`: [m/s] y velocity
  - `ax`: [m/s²] x acceleration
  - `ay`: [m/s²] y acceleration
  - `phi`: [rad] rotation angle (counterclockwise)
  - `ω`: [rad/s] angular velocity
  - `α`: [rad/s²] angular acceleration

# Connectors:

  - `frame`: 2-dim. Coordinate system

https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Parts/Body.mo
"""
@component function Body(; name, m, j, rx = 0, ry = 0, gy = -9.807)
    @named frame = Frame()
    pars = @parameters begin
        m = m
        j = j
        gy = gy
    end

    vars = @variables begin
        fx(t) = 0
        fy(t) = 0
        rx(t) = rx
        ry(t) = ry
        vx(t) = 0
        vy(t) = 0
        ax(t) = 0
        ay(t) = 0
        phi(t) = 0
        ω(t) = 0
        α(t) = 0
    end

    eqs = [
        # velocity is the time derivative of position
        rx ~ frame.x,
        ry ~ frame.y,
        vx ~ D(rx),
        vy ~ D(ry),
        phi ~ frame.phi,
        ω ~ D(phi),
        # acceleration is the time derivative of velocity
        ax ~ D(vx),
        ay ~ D(vy),
        α ~ D(ω),
        # newton's law
        fx ~ frame.fx,
        fy ~ frame.fy,
        ax ~ fx / m,
        ay ~ ifelse(gy !== nothing, fy / m + gy, fy / m),
        j * α ~ frame.j,
    ]

    return compose(ODESystem(eqs, t, vars, pars; name = name),
        frame)
end

"""
    FixedTranslation(; name, r::AbstractArray, l)
A fixed translation between two components (rigid rod)

# Parameters:

  - `rx`: [m] Fixed x-length of the rod resolved w.r.t to body frame_a at phi = 0
  - `ry`: [m] Fixed y-length of the rod resolved w.r.t to body frame_a at phi = 0

# Connectors:

      - `frame_a` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
      - `frame_b` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque

https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Parts/FixedTranslation.mo
"""
@mtkmodel FixedTranslation begin
    @extend frame_a, frame_b = partial_frames = PartialTwoFrames()

    @parameters begin
        rx = 0,
        [
            description = "Fixed x-length of the rod resolved w.r.t to body frame_a at phi = 0",
        ]
        ry = 0,
        [
            description = "Fixed y-length of the rod resolved w.r.t to body frame_a at phi = 0",
        ]
    end

    begin
        R = [cos(frame_a.phi) -sin(frame_a.phi);
            sin(frame_a.phi) cos(frame_a.phi)]
        r0 = R * [rx, ry]
    end

    @equations begin
        # rigidly connect positions
        frame_a.x + rx ~ frame_b.x
        frame_a.y + ry ~ frame_b.y
        frame_a.phi ~ frame_b.phi
        # balancing force including lever principle
        frame_a.fx + frame_b.fx ~ 0
        frame_a.fy + frame_b.fy ~ 0
        frame_a.j + frame_b.j + r0' * [frame_b.fy, -frame_b.fx] ~ 0
    end
end
