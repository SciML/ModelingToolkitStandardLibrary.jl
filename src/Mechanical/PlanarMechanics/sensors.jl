"""
    PartialAbsoluteSensor(;name)
Partial absolute sensor model for sensors defined by components
# Connectors:

  - `frame: 2-dim. Coordinate system
"""
@mtkmodel PartialAbsoluteSensor begin
    @components begin
        frame_a = Frame()
    end
    # TODO: assert the number of connections
    # https://github.com/dzimmer/PlanarMechanics/blob/443b007bcc1522bb172f13012e2d7a8ecc3f7a9b/PlanarMechanics/Sensors/Internal/PartialAbsoluteSensor.mo#L11
end

"""
    PartialRelativeSensor(;name)
Partial relative sensor model for sensors defined by components

# Connectors:

  - `frame_a`: Coordinate system a
  - `frame_b`: Coordinate system b
"""
@mtkmodel PartialRelativeSensor begin
    @components begin
        frame_a = Frame()
        frame_b = Frame()
    end
    # TODO: assert the number of connections
    # https://github.com/dzimmer/PlanarMechanics/blob/443b007bcc1522bb172f13012e2d7a8ecc3f7a9b/PlanarMechanics/Sensors/Internal/PartialRelativeSensor.mo#L12-L13
end

"""
    PartialAbsoluteBaseSensor(;name)
Partial absolute sensor models for sensors defined by equations (frame_resolve must be connected exactly once)
# Connectors:

  - `frame_a`: 2-dim. Coordinate system from which kinematic quantities are measured
  - `frame_resolve`: 2-dim. Coordinate system in which vector is optionally resolved
"""
@mtkmodel PartialAbsoluteBaseSensor begin
    @components begin
        frame_a = Frame()
        frame_resolve = FrameResolve()
    end

    @equations begin
        # TODO: assert the number of connections
        # https://github.com/dzimmer/PlanarMechanics/blob/443b007bcc1522bb172f13012e2d7a8ecc3f7a9b/PlanarMechanics/Sensors/Internal/PartialAbsoluteBaseSensor.mo#L20-L21
        frame_a.fx ~ 0
        frame_a.fy ~ 0
        frame_a.j ~ 0
        frame_resolve.fx ~ 0
        frame_resolve.fy ~ 0
        frame_resolve.j ~ 0
    end
end

"""
    PartialRelativeBaseSensor(;name)
Partial relative sensor models for sensors defined by equations (frame_resolve must be connected exactly once)

# Connectors:
  - `frame_a`: 
  - `frame_b`: 
  - `frame_resolve`: 
"""
@mtkmodel PartialRelativeBaseSensor begin
    @components begin
        frame_a = Frame()
        frame_b = Frame()
        frame_resolve = FrameResolve()
    end

    @equations begin
        # TODO: assert the number of connections
        # https://github.com/dzimmer/PlanarMechanics/blob/443b007bcc1522bb172f13012e2d7a8ecc3f7a9b/PlanarMechanics/Sensors/Internal/PartialRelativeBaseSensor.mo#L19-L21

        frame_a.fx ~ 0
        frame_a.fy ~ 0
        frame_a.j ~ 0
        frame_b.fx ~ 0
        frame_b.fy ~ 0
        frame_b.j ~ 0
        frame_resolve.fx ~ 0
        frame_resolve.fy ~ 0
        frame_resolve.j ~ 0
    end
end

"""
    BasicAbsolutePosition(;name, resolve_in_frame = :frame_a)
Measure absolute position and orientation (same as Sensors.AbsolutePosition, but frame_resolve is not conditional and must be connected).

# Connectors:

  - `x`: [m] x-position
  - `y`: [m] y-position
  - `phi`: [rad] rotation angle (counterclockwise)
  - `frame_a`: Coordinate system a
  - `frame_resolve`: 2-dim. Coordinate system in which vector is optionally resolved

# Parameters:
  
  - `resolve_in_frame`: Frame in which output x, y, phi r is resolved (1: :world, 2: :frame_a, 3: :frame_resolve)
"""
@component function BasicAbsolutePosition(; name, resolve_in_frame = :frame_a)
    @named x = RealOutput()
    @named y = RealOutput()
    @named phi = RealOutput()

    @named partial_abs_base_sensor = PartialAbsoluteBaseSensor()
    @unpack frame_a, frame_resolve = partial_abs_base_sensor

    if resolve_in_frame == :world
        r = [
            frame_a.x,
            frame_a.y,
            frame_a.phi]
    elseif resolve_in_frame == :frame_a
        rotation_matrix = [cos(frame_a.phi) -sin(frame_a.phi) 0;
                           sin(frame_a.phi) cos(frame_a.phi) 0;
                           0 0 1]
        r = transpose(rotation_matrix) * [frame_a.x, frame_a.y, frame_a.phi] -
            [0, 0, frame_a.phi]
    elseif resolve_in_frame == :frame_resolve
        rotation_matrix = [cos(frame_resolve.phi) -sin(frame_resolve.phi) 0;
                           sin(frame_resolve.phi) cos(frame_resolve.phi) 0;
                           0 0 1]
        r = transpose(rotation_matrix) * [frame_a.x, frame_a.y, frame_a.phi] -
            [0, 0, frame_resolve.phi]
    else
        error("resolve_in_frame must be one of :world, :frame_a, :frame_resolve")
    end

    eqs = [
        x.u ~ r[1],
        y.u ~ r[2],
        phi.u ~ r[3],
        frame_a.fx ~ 0,
        frame_a.fy ~ 0,
        frame_a.j ~ 0
    ]

    return compose(ODESystem(eqs, t, [], []; name = name),
        x, y, phi, frame_a, frame_resolve)
end

"""
    AbsolutePosition(;name, resolve_in_frame = :frame_a)
Measure absolute position and orientation of the origin of frame connector

# Connectors:

  - `x`: [m] x-position
  - `y`: [m] y-position
  - `phi`: [rad] rotation angle (counterclockwise)

# Parameters:

  - `resolve_in_frame`: Frame in which output x, y, phi is resolved (1: :world, 2: :frame_a, 3: :frame_resolve)
"""
@component function AbsolutePosition(; name, resolve_in_frame = :frame_a)
    @named pos = BasicAbsolutePosition(; resolve_in_frame)
    @named partial_abs_sensor = PartialAbsoluteSensor()
    @unpack frame_a, = partial_abs_sensor

    @named x = RealOutput()
    @named y = RealOutput()
    @named phi = RealOutput()

    systems = [pos, frame_a, x, y, phi]

    eqs = [
        x.u ~ pos.x.u,
        y.u ~ pos.y.u,
        phi.u ~ pos.phi.u,
        connect(pos.frame_a, frame_a)
    ]

    if resolve_in_frame == :frame_resolve
        @named fr = FrameResolve()
        push!(systems, fr)
        push!(eqs, connect(pos.frame_resolve, fr))
    end

    if resolve_in_frame != :frame_resolve
        @named zero_position = ZeroPosition()
        push!(systems, zero_position)
        push!(eqs, connect(zero_position.frame_resolve, pos.frame_resolve))
    end

    return compose(ODESystem(eqs, t, [], []; name = name), systems...)
end

"""
    BasicRelativePosition(; name, resolve_in_frame = :frame_a)
Measure relative position and orientation between the origins of two frame connectors

# Connectors:

 - `rel_x`: [m] Relative x-position
 - `rel_y`: [m] Relative y-position
 - `rel_phi`: [rad] Relative rotation angle (counterclockwise)
 - `frame_a`: Coordinate system a
 - `frame_b`: Coordinate system b
 - `frame_resolve`: 

# Parameters:
    
    - `resolve_in_frame`: Frame in which output x, y, phi is resolved (1: :world, 2: :frame_a, 3: frame_b 4: :frame_resolve)
"""
@component function BasicRelativePosition(; name, resolve_in_frame = :frame_a)
    @named rel_x = RealOutput()
    @named rel_y = RealOutput()
    @named rel_phi = RealOutput()

    @named partial_rel_pos = PartialRelativeBaseSensor()
    @unpack frame_a, frame_b, frame_resolve = partial_rel_pos

    if resolve_in_frame == :frame_a
        rotation_matrix = [cos(frame_a.phi) -sin(frame_a.phi) 0;
                           sin(frame_a.phi) cos(frame_a.phi) 0;
                           0 0 1]
        r = transpose(rotation_matrix) *
            [frame_b.x - frame_a.x, frame_b.y - frame_a.y, frame_b.phi - frame_a.phi]
    elseif resolve_in_frame == :frame_b
        rotation_matrix = [cos(frame_b.phi) -sin(frame_b.phi) 0;
                           sin(frame_b.phi) cos(frame_b.phi) 0;
                           0 0 1]
        r = transpose(rotation_matrix) *
            [frame_b.x - frame_a.x, frame_b.y - frame_a.y, frame_b.phi - frame_a.phi]
    elseif resolve_in_frame == :world
        r = [frame_b.x - frame_a.x, frame_b.y - frame_a.y, frame_b.phi - frame_a.phi]
    elseif resolve_in_frame == :frame_resolve
        rotation_matrix = [cos(frame_resolve.phi) -sin(frame_resolve.phi) 0;
                           sin(frame_resolve.phi) cos(frame_resolve.phi) 0;
                           0 0 1]
        r = transpose(rotation_matrix) *
            [frame_b.x - frame_a.x, frame_b.y - frame_a.y, frame_b.phi - frame_a.phi]
    else
        error("resolve_in_frame must be one of :world, :frame_a, :frame_resolve")
    end

    eqs = [
        rel_x.u ~ r[1],
        rel_y.u ~ r[2],
        rel_phi.u ~ r[3],
        frame_a.fx ~ 0,
        frame_a.fy ~ 0,
        frame_a.j ~ 0,
        frame_b.fx ~ 0,
        frame_b.fy ~ 0,
        frame_b.j ~ 0
    ]

    return compose(ODESystem(eqs, t, [], []; name = name),
        rel_x, rel_y, rel_phi, frame_a, frame_b, frame_resolve)
end

"""
    RelativePosition(; name, resolve_in_frame = :frame_a)
Measure relative position and orientation between the origins of two frame connectors

# Connectors:

    - `rel_x`: [m] Relative x-position
    - `re_y`: [m] Relative y-position
    - `rel_phi`: [rad] Relative rotation angle (counterclockwise)
    - `frame_a`: Coordinate system a
    - `frame_b`: Coordinate system b
# Parameters:

    - `resolve_in_frame`: Frame in which output x, y, phi is resolved (1: :world, 2: :frame_a, 3: frame_b 4: :frame_resolve)
"""
@component function RelativePosition(; name, resolve_in_frame = :frame_a)
    @named pos = BasicRelativePosition(; resolve_in_frame)
    @named partial_rel_pos = PartialRelativeSensor()
    @unpack frame_a, frame_b = partial_rel_pos

    @named rel_x = RealOutput()
    @named rel_y = RealOutput()
    @named rel_phi = RealOutput()

    systems = [pos, frame_a, frame_b, rel_x, rel_y, rel_phi]
    eqs = [
        pos.rel_x.u ~ rel_x.u,
        pos.rel_y.u ~ rel_y.u,
        pos.rel_phi.u ~ rel_phi.u,
        connect(pos.frame_a, frame_a),
        connect(pos.frame_b, frame_b)
    ]

    if resolve_in_frame == :frame_resolve
        @named fr = FrameResolve()
        push!(systems, fr)
        push!(eqs, connect(pos.frame_resolve, fr))
    end

    if resolve_in_frame != :frame_resolve
        @named zero_position = ZeroPosition()
        push!(systems, zero_position)
        push!(eqs, connect(zero_position.frame_resolve, pos.frame_resolve))
    end

    return compose(ODESystem(eqs, t, [], []; name = name),
        systems...)
end

@component function BasicTransformAbsoluteVector(;
        name,
        frame_in = :frame_a,
        frame_out = frame_in)
    @named x_in = RealInput()
    @named y_in = RealInput()
    @named phi_in = RealInput()

    @named x_out = RealOutput()
    @named y_out = RealOutput()
    @named phi_out = RealOutput()

    @named frame_a = Frame()
    @named frame_b = Frame()
    @named frame_resolve = FrameResolve()

    systems = [frame_a, frame_b, frame_resolve, x_in, y_in, phi_in, x_out, y_out, phi_out]
    eqs = [
        # TODO: assert the number of connections
        # https://github.com/dzimmer/PlanarMechanics/blob/443b007bcc1522bb172f13012e2d7a8ecc3f7a9b/PlanarMechanics/Sensors/Internal/BasicTransformAbsoluteVector.mo#L42-L43
        frame_a.fx ~ 0,
        frame_a.fy ~ 0,
        frame_a.j ~ 0,
        frame_resolve.fx ~ 0,
        frame_resolve.fy ~ 0,
        frame_resolve.j ~ 0
    ]

    r_temp = Vector{Float64}(undef, 3)
    R1 = Matrix{Float64}(undef, 3, 4)

    if frame_out == frame_in
        append!(eqs, [
            x_out.u ~ x_in.u,
            y_out.u ~ y_in.u,
            phi_out.u ~ phi_in.u
        ])
    else
        if frame_in == :world
            R1 = [1.0 0.0 0.0 0.0;
                  0.0 1.0 0.0 0.0;
                  0.0 0.0 1.0 0.0]
        elseif frame_in == :frame_a
            R1 = [cos(frame_a.phi) -sin(frame_a.phi) 0.0 0.0;
                  sin(frame_a.phi) cos(frame_a.phi) 0.0 0.0;
                  0.0 0.0 1.0 frame_a.phi]
        elseif frame_in == :frame_resolve
            R1 = [cos(frame_resolve.phi) -sin(frame_resolve.phi) 0.0 0.0;
                  sin(frame_resolve.phi) cos(frame_resolve.phi) 0.0 0.0;
                  0.0 0.0 1.0 frame_resolve.phi]
        else
            error("Wrong value for parameter frame_in")
        end

        r_temp = R1 * [x_in.u, y_in.u, phi_in.u, 1]

        if frame_out == :world
            append!(eqs, [
                x_out.u ~ r_temp[1],
                y_out.u ~ r_temp[2],
                phi_out.u ~ r_temp[3]
            ])
        elseif frame_out == :frame_a
            rotation_matrix = [cos(frame_a.phi) sin(frame_a.phi) 0.0;
                               -sin(frame_a.phi) cos(frame_a.phi) 0.0;
                               0.0 0.0 1.0]
            r = rotation_matrix * r_temp
            append!(eqs, [
                x_out.u ~ r[1],
                y_out.u ~ r[2],
                phi_out.u ~ r[3]
            ])
        elseif frame_out == :frame_resolve
            rotation_matrix = [cos(frame_resolve.phi) sin(frame_resolve.phi) 0.0;
                               -sin(frame_resolve.phi) cos(frame_resolve.phi) 0.0;
                               0.0 0.0 1.0]
            r = rotation_matrix * r_temp
            append!(eqs, [
                x_out.u ~ r[1],
                y_out.u ~ r[2],
                phi_out.u ~ r[3]
            ])
        else
            error("Wrong value for parameter frame_out")
        end
    end

    return compose(ODESystem(eqs, t, [], []; name = name), systems...)
end

@component function TransformAbsoluteVector(;
        name,
        frame_in = :frame_a,
        frame_out = frame_in)
    @named frame_a = Frame()

    @named x_in = RealInput()
    @named y_in = RealInput()
    @named phi_in = RealInput()

    @named x_out = RealOutput()
    @named y_out = RealOutput()
    @named phi_out = RealOutput()
    @named basic_transform_vector = BasicTransformAbsoluteVector(; frame_in, frame_out)

    systems = [frame_a, x_in, y_in, phi_in, x_out, y_out, phi_out, basic_transform_vector]

    eqs = [
        connect(basic_transform_vector.frame_a, frame_a),
        # out
        connect(basic_transform_vector.x_out, x_out),
        connect(basic_transform_vector.y_out, y_out),
        connect(basic_transform_vector.phi_out, phi_out),
        # in
        connect(basic_transform_vector.x_in, x_in),
        connect(basic_transform_vector.y_in, y_in),
        connect(basic_transform_vector.phi_in, phi_in)
    ]

    if frame_in == :frame_resolve || frame_out == :frame_resolve
        @named frame_resolve = FrameResolve()
        push!(systems, frame_resolve)
        push!(eqs, connect(basic_transform_vector.frame_resolve, frame_resolve))
    end

    if !(frame_in == :frame_resolve || frame_out == :frame_resolve)
        @named zero_pos = ZeroPosition()
        push!(systems, zero_pos)
        push!(eqs,
            connect(zero_pos.frame_resolve, basic_transform_vector.frame_resolve))
    end

    return compose(ODESystem(eqs, t, [], []; name = name), systems...)
end

@component function AbsoluteVelocity(; name, resolve_in_frame = :frame_a)
    @named partial_abs_sensor = PartialAbsoluteSensor()
    @unpack frame_a, = partial_abs_sensor

    @named v_x = RealOutput()
    @named v_y = RealOutput()
    @named ω = RealOutput()

    @named pos = BasicAbsolutePosition(; resolve_in_frame = :world)
    @named zero_pos = ZeroPosition()

    @named transform_absolute_vector = TransformAbsoluteVector(;
        frame_in = :world,
        frame_out = resolve_in_frame)

    systems = [frame_a, v_x, v_y, ω, pos, zero_pos, transform_absolute_vector]

    eqs = [
        # connect(position.r, der1.u),
        # connect(der1.y, transformAbsoluteVector.r_in)
        D(pos.x.u) ~ transform_absolute_vector.x_in.u,
        D(pos.y.u) ~ transform_absolute_vector.y_in.u,
        D(pos.phi.u) ~ transform_absolute_vector.phi_in.u,
        # connect(transformAbsoluteVector.r_out, v)
        transform_absolute_vector.x_out.u ~ v_x.u,
        transform_absolute_vector.y_out.u ~ v_y.u,
        transform_absolute_vector.phi_out.u ~ ω.u,
        connect(pos.frame_a, frame_a),
        connect(zero_pos.frame_resolve, pos.frame_resolve),
        connect(transform_absolute_vector.frame_a, frame_a)
    ]

    if resolve_in_frame == :frame_resolve
        @named frame_resolve = FrameResolve()
        push!(systems, frame_resolve)
        push!(eqs, connect(transform_absolute_vector.frame_resolve, frame_resolve))
    end

    if resolve_in_frame != :frame_resolve
        @named zero_pos1 = ZeroPosition()
        push!(systems, zero_pos1)
        push!(eqs,
            connect(transform_absolute_vector.zero_pos.frame_resolve,
                zero_pos1.frame_resolve))
    end

    return compose(ODESystem(eqs, t, [], []; name = name), systems...)
end

@component function BasicTransformRelativeVector(;
        name,
        frame_in = :frame_a,
        frame_out = frame_in)
    @named x_in = RealInput()
    @named y_in = RealInput()
    @named phi_in = RealInput()

    @named x_out = RealOutput()
    @named y_out = RealOutput()
    @named phi_out = RealOutput()

    @named partial_relative_base_sensor = PartialRelativeBaseSensor()
    @unpack frame_a, frame_b, frame_resolve = partial_relative_base_sensor

    systems = [
        x_in,
        y_in,
        phi_in,
        x_out,
        y_out,
        phi_out,
        frame_a,
        frame_b,
        frame_resolve
    ]

    eqs = Equation[]

    R1 = Matrix{Float64}(undef, 3, 3)

    if frame_out == frame_in
        append!(eqs, [
            x_out.u ~ x_in.u,
            y_out.u ~ y_in.u,
            phi_out.u ~ phi_in.u
        ])
    else
        if frame_in == :world
            R1 = [1.0 0.0 0.0;
                  0.0 1.0 0.0;
                  0.0 0.0 1.0]
        elseif frame_in == :frame_a
            R1 = [cos(frame_a.phi) -sin(frame_a.phi) 0.0;
                  sin(frame_a.phi) cos(frame_a.phi) 0.0;
                  0.0 0.0 1.0]
        elseif frame_in == :frame_b
            R1 = [cos(frame_b.phi) -sin(frame_b.phi) 0.0;
                  sin(frame_b.phi) cos(frame_b.phi) 0.0;
                  0.0 0.0 1.0]
        else
            R1 = [cos(frame_resolve.phi) -sin(frame_resolve.phi) 0.0;
                  sin(frame_resolve.phi) cos(frame_resolve.phi) 0.0;
                  0.0 0.0 1.0]
        end

        r_in = [x_in.u, y_in.u, phi_in.u]
        if frame_out == :world
            r_out = R1 * r_in
        elseif frame_out == :frame_a
            rotation_matrix = [cos(frame_a.phi) -sin(frame_a.phi) 0.0;
                               sin(frame_a.phi) cos(frame_a.phi) 0.0;
                               0.0 0.0 1.0]
            r_out = transpose(rotation_matrix) * (R1 * r_in)
        elseif frame_out == :frame_b
            rotation_matrix = [cos(frame_b.phi) -sin(frame_b.phi) 0.0;
                               sin(frame_b.phi) cos(frame_b.phi) 0.0;
                               0.0 0.0 1.0]
            r_out = transpose(rotation_matrix) * (R1 * r_in)
        else
            rotation_matrix = [cos(frame_resolve.phi) -sin(frame_resolve.phi) 0.0;
                               sin(frame_resolve.phi) cos(frame_resolve.phi) 0.0;
                               0.0 0.0 1.0]
            r_out = transpose(rotation_matrix) * (R1 * r_in)
        end
        append!(eqs, [
            x_out.u ~ r_out[1],
            y_out.u ~ r_out[2],
            phi_out.u ~ r_out[3]
        ])
    end

    return compose(ODESystem(eqs, t, [], []; name = name), systems...)
end

@component function TransformRelativeVector(;
        name,
        frame_in = :frame_a,
        frame_out = frame_in)
    @named x_in = RealInput()
    @named y_in = RealInput()
    @named phi_in = RealInput()

    @named x_out = RealOutput()
    @named y_out = RealOutput()
    @named phi_out = RealOutput()

    @named basic_transformb_vector = BasicTransformRelativeVector(; frame_in, frame_out)
    @named partial_relative_sensor = PartialRelativeSensor()
    @unpack frame_a, frame_b = partial_relative_sensor

    systems = [
        x_in,
        y_in,
        phi_in,
        x_out,
        y_out,
        phi_out,
        basic_transformb_vector,
        frame_a,
        frame_b
    ]

    eqs = [
        connect(basic_transformb_vector.frame_a, frame_a),
        connect(basic_transformb_vector.frame_b, frame_b),
        # out
        connect(basic_transformb_vector.x_out, x_out),
        connect(basic_transformb_vector.y_out, y_out),
        connect(basic_transformb_vector.phi_out, phi_out),
        # in
        connect(basic_transformb_vector.x_in, x_in),
        connect(basic_transformb_vector.y_in, y_in),
        connect(basic_transformb_vector.phi_in, phi_in)
    ]

    if frame_in == :frame_resolve || frame_out == :frame_resolve
        @named frame_resolve = FrameResolve()
        push!(systems, frame_resolve)
        push!(eqs, connect(basic_transformb_vector.frame_resolve, frame_resolve))
    end

    if !(frame_in == :frame_resolve || frame_out == :frame_resolve)
        @named zero_pos = ZeroPosition()
        push!(systems, zero_pos)
    end

    return compose(ODESystem(eqs, t, [], []; name = name), systems...)
end

@component function RelativeVelocity(; name, resolve_in_frame = :frame_a)
    @named partial_rel_sensor = PartialRelativeSensor()
    @unpack frame_a, frame_b = partial_rel_sensor

    @named rel_v_x = RealOutput()
    @named rel_v_y = RealOutput()
    @named rel_ω = RealOutput()

    @named rel_pos = RelativePosition(; resolve_in_frame = :frame_a)
    @named transform_relative_vector = TransformRelativeVector(;
        frame_in = :frame_a,
        frame_out = resolve_in_frame)

    systems = [
        frame_a,
        frame_b,
        rel_v_x,
        rel_v_y,
        rel_ω,
        rel_pos,
        transform_relative_vector
    ]
    eqs = [
        # connect(relativePosition.r_rel, der_r_rel.u) 
        # connect(der_r_rel.y, transformRelativeVector.r_in)
        D(rel_pos.rel_x.u) ~ transform_relative_vector.x_in.u,
        D(rel_pos.rel_y.u) ~ transform_relative_vector.y_in.u,
        D(rel_pos.rel_phi.u) ~ transform_relative_vector.phi_in.u,
        #  connect(transformRelativeVector.r_out, v_rel)
        transform_relative_vector.x_out.u ~ rel_v_x.u,
        transform_relative_vector.y_out.u ~ rel_v_y.u,
        transform_relative_vector.phi_out.u ~ rel_ω.u,
        connect(rel_pos.frame_a, frame_a),
        connect(rel_pos.frame_b, frame_b),
        connect(transform_relative_vector.frame_a, frame_a),
        connect(transform_relative_vector.frame_b, frame_b)
    ]

    if resolve_in_frame == :frame_resolve
        @named frame_resolve = FrameResolve()
        push!(systems, frame_resolve)
        push!(eqs, connect(transform_relative_vector.frame_resolve, frame_resolve))
    end

    if resolve_in_frame != :frame_resolve
        @named zero_pos = ZeroPosition()
        push!(systems, zero_pos)
        push!(eqs,
            connect(transform_relative_vector.zero_pos.frame_resolve,
                zero_pos.frame_resolve))
    end

    return compose(ODESystem(eqs, t, [], []; name = name), systems...)
end

@component function AbsoluteAcceleration(; name, resolve_in_frame = :frame_a)
    @named partial_abs_sensor = PartialAbsoluteSensor()
    @unpack frame_a, = partial_abs_sensor

    @named a_x = RealOutput()
    @named a_y = RealOutput()
    @named α = RealOutput()

    @named pos = BasicAbsolutePosition(; resolve_in_frame = :world)
    @named zero_pos = ZeroPosition()

    @named transform_absolute_vector = TransformAbsoluteVector(;
        frame_in = :world,
        frame_out = resolve_in_frame)

    systems = [frame_a, a_x, a_y, α, pos, zero_pos, transform_absolute_vector]

    eqs = [
        # connect(position.r, der1.u)
        # connect(der1.y, der2.u)
        # connect(der2.y, transformAbsoluteVector.r_in)
        D(D(pos.x.u)) ~ transform_absolute_vector.x_in.u,
        D(D(pos.y.u)) ~ transform_absolute_vector.y_in.u,
        D(D(pos.phi.u)) ~ transform_absolute_vector.phi_in.u,
        # connect(transformAbsoluteVector.r_out, a) 
        transform_absolute_vector.x_out.u ~ a_x.u,
        transform_absolute_vector.y_out.u ~ a_y.u,
        transform_absolute_vector.phi_out.u ~ α.u,
        connect(pos.frame_a, frame_a),
        connect(zero_pos.frame_resolve, pos.frame_resolve),
        connect(transform_absolute_vector.frame_a, frame_a)
    ]

    if resolve_in_frame == :frame_resolve
        @named frame_resolve = FrameResolve()
        push!(systems, frame_resolve)
        push!(eqs, connect(transform_absolute_vector.frame_resolve, frame_resolve))
    end

    if resolve_in_frame != :frame_resolve
        @named zero_pos1 = ZeroPosition()
        push!(systems, zero_pos1)
        push!(eqs,
            connect(transform_absolute_vector.zero_pos.frame_resolve,
                zero_pos1.frame_resolve))
    end

    return compose(ODESystem(eqs, t, [], []; name = name), systems...)
end

@component function RelativeAcceleration(; name, resolve_in_frame = :frame_a)
    @named partial_rel_sensor = PartialRelativeSensor()
    @unpack frame_a, frame_b = partial_rel_sensor

    @named rel_a_x = RealOutput()
    @named rel_a_y = RealOutput()
    @named rel_α = RealOutput()

    @named rel_pos = RelativePosition(; resolve_in_frame = :frame_a)
    @named transform_relative_vector = TransformRelativeVector(;
        frame_in = :frame_a,
        frame_out = resolve_in_frame)

    systems = [
        frame_a,
        frame_b,
        rel_a_x,
        rel_a_y,
        rel_α,
        rel_pos,
        transform_relative_vector
    ]

    eqs = [
        # connect(der_r_rel.y, transformRelativeVector.r_in)
        # connect(transformRelativeVector.r_out, der_r_rel1.u)
        # connect(relativePosition.r_rel, der_r_rel.u)
        D(D(rel_pos.rel_x.u)) ~ transform_relative_vector.x_in.u,
        D(D(rel_pos.rel_y.u)) ~ transform_relative_vector.y_in.u,
        D(D(rel_pos.rel_phi.u)) ~ transform_relative_vector.phi_in.u,
        # connect(der_r_rel1.y, a_rel),
        transform_relative_vector.x_out.u ~ rel_a_x.u,
        transform_relative_vector.y_out.u ~ rel_a_y.u,
        transform_relative_vector.phi_out.u ~ rel_α.u,
        connect(rel_pos.frame_a, frame_a),
        connect(rel_pos.frame_b, frame_b),
        connect(transform_relative_vector.frame_a, frame_a),
        connect(transform_relative_vector.frame_b, frame_b)
    ]

    if resolve_in_frame == :frame_resolve
        @named frame_resolve = FrameResolve()
        push!(systems, frame_resolve)
        push!(eqs, connect(transform_relative_vector.frame_resolve, frame_resolve))
    end

    if resolve_in_frame != :frame_resolve
        @named zero_pos = ZeroPosition()
        push!(systems, zero_pos)
        push!(eqs,
            connect(transform_relative_vector.zero_pos.frame_resolve,
                zero_pos.frame_resolve))
    end

    return compose(ODESystem(eqs, t, [], []; name = name), systems...)
end

function connect_sensor(component_frame, sensor_frame)
    # TODO: make this an override of the `connect` method
    return [
        component_frame.x ~ sensor_frame.x,
        component_frame.y ~ sensor_frame.y,
        component_frame.phi ~ sensor_frame.phi
    ]
end