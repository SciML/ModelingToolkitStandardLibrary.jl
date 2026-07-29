"""
    MechanicalPort(; name)

One-dimensional translational mechanical port.

# States
- `v(t)`: [`m/s`] Velocity of the node.
- `f(t)`: [`N`] Force entering the node.
"""
@connector function MechanicalPort(; name, v = nothing, f = nothing)
    vars = @variables begin
        v(t) = v
        f(t) = f, [connect = Flow]
    end
    return System(Equation[], t, vars, []; name)
end
