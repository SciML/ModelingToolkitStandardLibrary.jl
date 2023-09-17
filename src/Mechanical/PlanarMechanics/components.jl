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
