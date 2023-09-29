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
        x ~ r[1],
        y ~ r[2],
        phi ~ r[3],
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
    @named pos = BasicAbsolutePosition()
    @named partial_abs_sensor = PartialAbsoluteSensor()
    @unpack frame_a, = partial_abs_sensor

    @named x = RealOutput()
    @named y = RealOutput()
    @named phi = RealOutput()

    systems = [pos, frame_a, x, y, phi]

    eqs = [
        pos.x ~ x,
        pos.y ~ y,
        pos.phi ~ phi,
        connect(pos.frame_a, frame_a),
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

"""
    BasicRelativePosition(; name, resolve_in_frame = :frame_a)
Measure relative position and orientation between the origins of two frame connectors

# Connectors:

 - `x`: [m] x-position
 - `y`: [m] y-position
 - `phi`: [rad] rotation angle (counterclockwise)
 - `frame_a`: Coordinate system a
 - `frame_b`: Coordinate system b
 - `frame_resolve`: 

# Parameters:
    
    - `resolve_in_frame`: Frame in which output x, y, phi is resolved (1: :world, 2: :frame_a, 3: frame_b 4: :frame_resolve)
"""
@component function BasicRelativePosition(; name, resolve_in_frame = :frame_a)
    @named x = RealOutput()
    @named y = RealOutput()
    @named phi = RealOutput()

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
        x ~ r[1],
        y ~ r[2],
        phi ~ r[3],
    ]

    return compose(ODESystem(eqs, t, [], []; name = name),
        x, y, phi, frame_a, frame_b, frame_resolve)
end

"""
    RelativePosition(; name, resolve_in_frame = :frame_a)
Measure relative position and orientation between the origins of two frame connectors

# Connectors:

    - `x`: [m] x-position
    - `y`: [m] y-position
    - `phi`: [rad] rotation angle (counterclockwise)
    - `frame_a`: Coordinate system a
    - `frame_b`: Coordinate system b
# Parameters:

    - `resolve_in_frame`: Frame in which output x, y, phi is resolved (1: :world, 2: :frame_a, 3: frame_b 4: :frame_resolve)
"""
@component function RelativePosition(; name, resolve_in_frame = :frame_a)
    @named pos = BasicRelativePosition(; resolve_in_frame)
    @named partial_rel_pos = PartialRelativeSensor()
    @unpack frame_a, frame_b = partial_rel_pos

    @named x = RealOutput()
    @named y = RealOutput()
    @named phi = RealOutput()

    systems = [pos, frame_a, frame_b, x, y, phi]
    eqs = [
        pos.x ~ x,
        pos.y ~ y,
        pos.phi ~ phi,
        connect(pos.frame_a, frame_a),
        connect(pos.frame_b, frame_b),
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
