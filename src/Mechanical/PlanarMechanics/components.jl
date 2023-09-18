"""
    Fixed(; name, r = (0.0, 0.0), phi = 0.0)

Frame fixed in the planar world frame at a given position and orientation

# Parameters:

  - `r`: [m, m] Fixed absolute x,y-position, resolved in planarWorld frame
  - `phi`: [rad] Fixed angle

# Connectors:

  - `frame: 2-dim. Coordinate system
"""
@mtkmodel Fixed begin
    @parameters begin
        r, [description = "Fixed absolute x,y-position, resolved in planarWorld frame"]
        phi, [description = "Fixed angle"]
    end

    @components begin
        frame = Frame(; r = (0.0, 0.0), phi = 0.0)
    end

    @equations begin
        frame.x, frame.y ~ r
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
  - `gy`: [m/s²] (optional) gravity field acting on the mass in the y-direction, positive value acts in the positive direction

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
"""
@component function Body(; name, m, j, r = nothing, gy = nothing)
    @named frame = Frame()
    pars = @parameters begin
        m = m
        j = j
        gy = gy
    end

    vars = @variables begin
        fx(t) = 0
        fy(t) = 0
        rx(t) = 0
        ry(t) = 0
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
        vx ~ D(r_x),
        vy ~ D(r_y),
        phi ~ frame.phi,
        ω ~ D(phi),
        # acceleration is the time derivative of velocity
        ax ~ D(v_x),
        ay ~ D(v_y),
        α ~ D(ω),
        # newton's law
        fx ~ frame.fx,
        fy ~ frame.fy,
        ax ~ fx / m,
        ay ~ ifelse(gy !== nothing, fy / m + gy, fy / m),
        j * α ~ frame.j,
    ]

    if r !== nothing
        rx = r[1]
        ry = r[2]
    end

    return compose(ODESystem(eqs, t, vars, pars; name = name),
        frame)
end
