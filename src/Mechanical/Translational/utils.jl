@connector function MechanicalPort(; name, f_int = 0, v_int = 0)
    pars = @parameters begin
        f_int = f_int
        v_int = v_int
    end
    vars = @variables begin
        v(t) = v_int
        f(t), [connect = Flow]
    end
    ODESystem(Equation[], t, vars, pars; name, defaults = [f => f_int])
end
Base.@doc """
    MechanicalPort(;name)

1-dim. rotational flange of a shaft.

# States:
- `v`: [m/s] velocity of the node
- `f`: [N] force entering the node
""" MechanicalPort
