@connector function MechanicalPort(; name, v = nothing, f = nothing)
    vars = @variables begin
        v(t) = v
        f(t) = f, [connect = Flow]
    end
    return System(Equation[], t, vars, []; name)
end
Base.@doc """
    MechanicalPort(;name)

1-dim. rotational flange of a shaft.

# States:
- `v`: [m/s] velocity of the node
- `f`: [N] force entering the node
""" MechanicalPort
