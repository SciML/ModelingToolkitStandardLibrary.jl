
@enum ResolveInFrameA begin
    world
    frame_a
    frame_resolve
end

Base.@doc """
    Enumeration to define the frame in which an absolute vector is resolved (world, frame_a, frame_resolve)
    Values:
        - `world`: Resolve in world frame
        - `frame_a`: Resolve in frame_a"
        - `frame_resolve`: Resolve in frame_resolve (frame_resolve must be connected)
""" ResolveInFrameA