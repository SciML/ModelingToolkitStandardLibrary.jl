@connector function RigidBody2DPort(; name)
    pars = []
    vars = @variables begin
        dx(t) = 0
        f_x(t), [connect = Flow]
        dy(t) = 0
        f_y(t), [connect = Flow]
        dA(t) = 0
        T_z(t), [connect = Flow]
    end
    ODESystem(Equation[], t, vars, pars, name = name, defaults = [f_x => 0, f_y => 0, T_z=>0])
end