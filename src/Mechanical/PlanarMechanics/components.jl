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
@component function Body(; name, m, I, rx = 0, ry = 0, gy = -9.807)
    @named frame = Frame()
    pars = @parameters begin
        m = m
        I = I
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
        I * α ~ frame.j,
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
        frame_a.x + r0[1] ~ frame_b.x
        frame_a.y + r0[2] ~ frame_b.y
        frame_a.phi ~ frame_b.phi
        # balancing force including lever principle
        frame_a.fx + frame_b.fx ~ 0
        frame_a.fy + frame_b.fy ~ 0
        frame_a.j + frame_b.j + r0' * [frame_b.fy, -frame_b.fx] ~ 0
    end
end

"""
    Spring(; name, c_x = 1, c_y = 1, c_phi = 1e5, s_relx0 = 0, s_rely0 = 0, phi_rel0 = 0, s_small = 1.e-10)

Linear 2D translational spring

# Parameters:

    - `c_x`: [N/m] Spring constant in x dir
    - `c_y`: [N/m] Spring constant in y dir
    - `c_phi`: [N.m/rad] Spring constant in phi dir
    - `s_relx0`: [m] Unstretched spring length
    - `s_rely0`: [m] Unstretched spring length
    - `phi_rel0`: [rad] Unstretched spring angle
    - `s_small`: [m] Prevent zero-division if distance between frame_a and frame_b is zero


# Connectors:

    - `frame_a` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
    - `frame_b` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque

https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Parts/Spring.mo
"""
@mtkmodel Spring begin
    @extend frame_a, frame_b = partial_frames = PartialTwoFrames()

    @parameters begin
        c_x = 1, [description = "Spring constant in x dir"]
        c_y = 1, [description = "Spring constant in y dir"]
        c_phi = 1.0e5, [description = "Spring constant"]
        s_relx0 = 0, [description = "Unstretched spring length"]
        s_rely0 = 0, [description = "Unstretched spring length"]
        phi_rel0 = 0, [description = "Unstretched spring angle"]
        s_small = 1.e-10,
        [
            description = "Prevent zero-division if distance between frame_a and frame_b is zero",
        ]
    end

    @variables begin
        s_relx(t) = 0
        s_rely(t) = 0
        phi_rel(t) = 0
        f_x(t)
        f_y(t)
    end

    begin
        r_rel_0 = [s_relx, s_rely, 0]
        l = sqrt(r_rel_0' * r_rel_0)
        e_rel_0 = r_rel_0 / max(l, s_small)
    end

    @equations begin
        phi_rel ~ frame_b.phi - frame_a.phi
        frame_a.j ~ 0
        frame_b.j ~ 0
        s_relx ~ frame_b.x - frame_a.x
        s_rely ~ frame_b.y - frame_a.y
        f_x ~ c_x * (s_relx - s_relx0)
        f_y ~ c_y * (s_rely - s_rely0)
        frame_a.fx ~ -f_x
        frame_b.fx ~ f_x
        frame_a.fy ~ -f_y
        frame_b.fy ~ f_y
    end
end

"""
    Damper(; name, d = 1, s_smal = 1.e-10)

Linear (velocity dependent) damper

# Parameters:

    - `d`: [N.s/m] Damoing constant 
    - `s_small`: [m] Prevent zero-division if distance between frame_a and frame_b is zero


# Connectors:

    - `frame_a` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
    - `frame_b` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque


https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Parts/Damper.mo
"""
@mtkmodel Damper begin
    @extend frame_a, frame_b = partial_frames = PartialTwoFrames()

    @parameters begin
        d = 1, [description = "damping constant"]
        s_small = 1.e-10,
        [
            description = "Prevent zero-division if distance between frame_a and frame_b is zero",
        ]
    end

    @variables begin
        r0x(t) = 0
        r0y(t) = 0
        d0x(t) = 0
        d0y(t) = 0
        vx(t) = 0
        vy(t) = 0
        v(t)
        f(t)
    end

    begin
        r0 = [r0x, r0y]
        l = sqrt(r0' * r0)
    end

    @equations begin
        frame_a.x + r0x ~ frame_b.x
        frame_a.y + r0y ~ frame_b.y
        D(frame_a.x) + vx ~ D(frame_b.x)
        D(frame_a.y) + vy ~ D(frame_b.y)
        v ~ [vx, vy]' * [d0x, d0y]
        f ~ -d * v
        d0x ~ ifelse(l < s_small, r0[1], r0[1] / l)
        d0y ~ ifelse(l < s_small, r0[2], r0[2] / l)
        frame_a.fx ~ d0x * f
        frame_a.fy ~ d0y * f
        frame_a.j ~ 0
        frame_a.fx + frame_b.fx ~ 0
        frame_a.fy + frame_b.fy ~ 0
        frame_a.j + frame_b.j ~ 0

        # lossPower ~ -f * v
    end
end

"""
    SpringDamper(; name, c_x = 1, c_y = 1, c_phi = 1e5, d_x = 1, d_y = 1, d_phi = 1, s_relx0 = 0, s_rely0 = 0, phi_rel0 = 0, s_small = 1.e-10)

Linear 2D translational spring damper model

# Parameters:

    - `c_x`: [N/m] Spring constant in x dir
    - `c_y`: [N/m] Spring constant in y dir
    - `c_phi`: [N.m/rad] Spring constant in phi dir
    - `d_x`: [N.s/m] Damping constant in x dir
    - `d_y`: [N.s/m] Damping constant in y dir
    - `d_phi`: [N.m.s/rad] Damping constant in phi dir
    - `s_relx0`: [m] Unstretched spring length
    - `s_rely0`: [m] Unstretched spring length
    - `phi_rel0`: [rad] Unstretched spring angle
    - `s_small`: [m] Prevent zero-division if distance between frame_a and frame_b is zero


# Connectors:

    - `frame_a` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
    - `frame_b` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque

https://github.com/dzimmer/PlanarMechanics/blob/743462f58858a808202be93b708391461cbe2523/PlanarMechanics/Parts/SpringDamper.mo
"""
@mtkmodel SpringDamper begin
    @extend frame_a, frame_b = partial_frames = PartialTwoFrames()

    @parameters begin
        c_x = 1, [description = "Spring constant in x dir"]
        c_y = 1, [description = "Spring constant in y dir"]
        c_phi = 1.0e5, [description = "Spring constant in phi dir"]
        d_x = 1, [description = "Damping constant in x dir"]
        d_y = 1, [description = "Damping constant in y dir"]
        d_phi = 1, [description = "Damping constant in phi dir"]
        s_relx0 = 0, [description = "Unstretched spring length"]
        s_rely0 = 0, [description = "Unstretched spring length"]
        phi_rel0 = 0, [description = "Unstretched spring angle"]
        s_small = 1.e-10,
        [
            description = "Prevent zero-division if distance between frame_a and frame_b is zero",
        ]
    end

    @variables begin
        v_relx(t)
        v_rely(t)
        ω_rel(t) = 0
        s_relx(t)
        s_rely(t)
        phi_rel(t) = 0
        f_x(t)
        f_y(t)
        j(t)
    end

    begin
        r_rel_0 = [s_relx, s_rely, 0]
        l = sqrt(r_rel_0' * r_rel_0)
        e_rel_0 = r_rel_0 / max(l, s_small)
    end

    @equations begin
        s_relx ~ frame_b.x - frame_a.x
        s_rely ~ frame_b.y - frame_a.y
        phi_rel ~ frame_b.phi - frame_a.phi
        v_relx ~ D(s_relx)
        v_rely ~ D(s_rely)
        ω_rel ~ D(phi_rel)

        j ~ c_phi * (phi_rel - phi_rel0) + d_phi * ω_rel
        frame_a.j ~ -j
        frame_b.j ~ j
        f_x ~ c_x * (s_relx - s_relx0) + d_x * v_relx
        f_y ~ c_y * (s_rely - s_rely0) + d_y * v_rely
        frame_a.fx ~ -f_x
        frame_b.fx ~ f_x
        frame_a.fy ~ -f_y
        frame_b.fy ~ f_y

        # lossPower ~ d_x * v_relx * v_relx + d_y * v_rely * v_rely
    end
end
