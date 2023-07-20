@connector MechanicalPort begin
    v(t) = 0.0
    f(t) = 0.0, [connect = Flow]
end
Base.@doc """
    MechanicalPort(;name)

1-dim. rotational flange of a shaft.

# States:
- `v`: [m/s] velocity of the node
- `f`: [N] force entering the node
""" MechanicalPort
