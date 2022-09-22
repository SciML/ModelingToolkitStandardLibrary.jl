@connector function Port(; name, v0 = 0.0)
    pars = @parameters v0 = v0
    vars = @variables begin
        v(t)
        f(t), [connect = Flow]
    end
    ODESystem(Equation[], t, vars, pars, name = name, defaults = Dict(v => v0, f => 0))
end
Base.@doc """
    Port(;name, v0=0.0, f0=0.0)

1-dim. rotational flange of a shaft.

# States:
- `v`: [m/s] velocity of the node
- `f`: [N] force entering the node
""" Port
