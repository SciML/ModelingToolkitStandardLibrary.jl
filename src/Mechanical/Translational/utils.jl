@connector MechanicalPort begin
    v(t) = 0.0, [description = "Velocity of the node", unit = u"m/s"]
    f(t) = 0.0, [connect = Flow, description = "Force entering the node", unit = u"N"]
end
Base.@doc """
    MechanicalPort(;name)

1-dim. rotational flange of a shaft.

# States:
- `v`: [m/s] velocity of the node
- `f`: [N] force entering the node
""" MechanicalPort
