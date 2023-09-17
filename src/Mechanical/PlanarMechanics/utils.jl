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