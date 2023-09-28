@connector Frame begin
    x(t), [description = "x position"]
    y(t), [description = "y position"]
    phi(t), [description = "rotation angle (counterclockwise)"]
    fx(t), [connect = Flow, description = "force in the x direction"]
    fy(t), [connect = Flow, description = "force in the y direction"]
    j(t), [connect = Flow, description = "torque (clockwise)"]
end
Base.@doc """
    Frame(;name)

Coordinate system (2-dim.) fixed to the component with one cut-force and cut-torque

# States:
    - `x`: [m] x position
    - `y`: [m] y position
    - `phi`: [rad] rotation angle (counterclockwise)
    - `fx`: [N] force in the x direction
    - `fy`: [N] force in the y direction
    - `j`: [N.m] torque (clockwise)
""" Frame

# extends Frame with just styling
# https://github.com/dzimmer/PlanarMechanics/blob/master/PlanarMechanics/Interfaces/Frame_resolve.mo
FrameResolve = Frame

@mtkmodel PartialTwoFrames begin
    @components begin
        frame_a = Frame()
        frame_b = Frame()
    end
end

Base.@doc """
    PartialTwoFrames(;name)
Partial model with two frames

# Connectors:
    - `frame_a` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
    - `frame_b` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
""" PartialTwoFrames
