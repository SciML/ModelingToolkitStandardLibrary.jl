@connector function Port(; name, v0=0.0, f0=0.0)
    sts = @variables begin
        v(t)
        f(t), [connect = Flow]
    end
    ODESystem(Equation[], t, sts, [], name = name, defaults = Dict(v => v0, f => f0))
end
Base.@doc """
    Port(;name, v_int=0.0, f_int=0.0)

1-dim. rotational flange of a shaft.

# States:
- `v`: [m/s] velocity of the node
- `f`: [N] force entering the node
""" Port
