@connector function Port(; name)
    pars = []
    vars = @variables begin
        v(t)
        f(t), [connect = Flow]
    end
    ODESystem(Equation[], t, vars, pars, name = name, defaults = [f => 0])
end
Base.@doc """
    Port(;name)

1-dim. rotational flange of a shaft.

# States:
- `v`: [m/s] velocity of the node
- `f`: [N] force entering the node
""" Port
