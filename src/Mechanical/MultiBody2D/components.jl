function Link(; name, m, l, I, g, x1_0 = 0, y1_0 = 0)
    pars = @parameters begin
        m = m
        l = l
        I = I
        g = g
        x1_0 = x1_0
        y1_0 = y1_0
    end

    vars = @variables begin
        (A(t) = 0), [state_priority = 10]
        (dA(t) = 0), [state_priority = 10]
        (ddA(t) = 0), [state_priority = 10]

        fx1(t) = 0
        fy1(t) = 0

        fx2(t) = 0
        fy2(t) = 0

        x1(t) = x1_0
        dx1(t) = 0

        y1(t) = y1_0
        dy1(t) = 0

        x2(t) = l + x1_0
        dx2(t) = 0

        y2(t) = 0
        dy2(t) = 0

        x_cm(t) = l / 2 + x1_0
        dx_cm(t) = 0
        ddx_cm(t) = 0

        y_cm(t) = 0
        dy_cm(t) = 0
        ddy_cm(t) = 0
    end

    @named TX1 = MechanicalPort()
    @named TY1 = MechanicalPort()

    @named TX2 = MechanicalPort()
    @named TY2 = MechanicalPort()

    eqs = [D(A) ~ dA
           D(dA) ~ ddA
           D(x1) ~ dx1
           D(y1) ~ dy1
           D(x2) ~ dx2
           D(y2) ~ dy2
           D(x_cm) ~ dx_cm
           D(dx_cm) ~ ddx_cm
           D(y_cm) ~ dy_cm
           D(dy_cm) ~ ddy_cm

           # x forces
           m * ddx_cm ~ fx1 + fx2

           # y forces
           m * ddy_cm ~ m * g + fy1 + fy2

           # torques
           I * ddA ~ -fy1 * (x2 - x1) / 2 + fy2 * (x2 - x1) / 2 + fx1 * (y2 - y1) / 2 -
                     fx2 * (y2 - y1) / 2

           # geometry
           x2 ~ l * cos(A) + x1
           y2 ~ l * sin(A) + y1
           x_cm ~ l * cos(A) / 2 + x1
           y_cm ~ l * sin(A) / 2 + y1
           TX1.f ~ fx1
           TX1.v ~ dx1
           TY1.f ~ fy1
           TY1.v ~ dy1
           TX2.f ~ fx2
           TX2.v ~ dx2
           TY2.f ~ fy2
           TY2.v ~ dy2]

    return ODESystem(eqs, t, vars, pars; name = name, systems = [TX1, TY1, TX2, TY2],
                     defaults = [TX1.v => 0, TY1.v => 0, TX2.v => 0, TY2.v => 0])
end
