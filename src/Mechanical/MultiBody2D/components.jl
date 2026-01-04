@component function Link(;
        name,
        m = nothing, l = nothing, I = nothing, g = nothing,
        x1_0 = 0.0, y1_0 = 0.0,
        A = nothing, dA = nothing, ddA = nothing,
        fx1 = nothing, fy1 = nothing,
        fx2 = nothing, fy2 = nothing,
        x1 = nothing, dx1 = nothing,
        y1 = nothing, dy1 = nothing,
        x2 = nothing, dx2 = nothing,
        y2 = nothing, dy2 = nothing,
        x_cm = nothing, dx_cm = nothing, ddx_cm = nothing,
        y_cm = nothing, dy_cm = nothing, ddy_cm = nothing
    )
    pars = @parameters begin
        m = m
        l = l
        I = I
        g = g
        x1_0 = x1_0
        y1_0 = y1_0
    end

    @named TX1 = Flange()
    @named TY1 = Flange()
    @named TX2 = Flange()
    @named TY2 = Flange()
    systems = [TX1, TY1, TX2, TY2]

    # Compute defaults using parameter symbols if not provided
    x1_val = x1 === nothing ? x1_0 : x1
    y1_val = y1 === nothing ? y1_0 : y1
    x2_val = x2 === nothing ? l + x1_0 : x2
    x_cm_val = x_cm === nothing ? l / 2 + x1_0 : x_cm

    vars = @variables begin
        A(t) = A, [state_priority = 10]
        dA(t) = dA, [state_priority = 10]
        ddA(t) = ddA, [state_priority = 10]

        fx1(t) = fx1
        fy1(t) = fy1

        fx2(t) = fx2
        fy2(t) = fy2

        x1(t) = x1_val
        dx1(t) = dx1

        y1(t) = y1_val
        dy1(t) = dy1

        x2(t) = x2_val
        dx2(t) = dx2

        y2(t) = y2
        dy2(t) = dy2

        x_cm(t) = x_cm_val
        dx_cm(t) = dx_cm
        ddx_cm(t) = ddx_cm

        y_cm(t) = y_cm
        dy_cm(t) = dy_cm
        ddy_cm(t) = ddy_cm
    end

    equations = Equation[
        D(A) ~ dA,
        D(dA) ~ ddA,
        D(x1) ~ dx1,
        D(y1) ~ dy1,
        D(x2) ~ dx2,
        D(y2) ~ dy2,
        D(x_cm) ~ dx_cm,
        D(dx_cm) ~ ddx_cm,
        D(y_cm) ~ dy_cm,
        D(dy_cm) ~ ddy_cm,

        # x forces
        m * ddx_cm ~ fx1 + fx2,

        # y forces
        m * ddy_cm ~ m * g + fy1 + fy2,

        # torques
        I * ddA ~
            -fy1 * (x2 - x1) / 2 + fy2 * (x2 - x1) / 2 + fx1 * (y2 - y1) / 2 -
            fx2 * (y2 - y1) / 2,

        # geometry
        x2 ~ l * cos(A) + x1,
        y2 ~ l * sin(A) + y1,
        x_cm ~ l * cos(A) / 2 + x1,
        y_cm ~ l * sin(A) / 2 + y1,
        TX1.f ~ fx1,
        TX1.s ~ x1,
        TY1.f ~ fy1,
        TY1.s ~ y1,
        TX2.f ~ fx2,
        TX2.s ~ x2,
        TY2.f ~ fy2,
        TY2.s ~ y2,
    ]

    return System(equations, t, vars, pars; name, systems)
end
