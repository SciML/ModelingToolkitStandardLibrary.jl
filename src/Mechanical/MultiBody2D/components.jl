function Link(; name, m, l, I, g)
    pars = @parameters begin
        m = m
        l = l
        I = I
        g = g
    end
    vars = @variables begin
        A(t) = 0

        x_1(t) = 0
        y_1(t) = 0

        x_2(t) = 0
        y_2(t) = 0

        fx1(t) = 0
        fy1(t) = 0

        fx2(t) = 0
        fy2(t) = 0
    end

    Port = Flange
    @named TX1 = Port()
    @named TY1 = Port()

    @named TX2 = Port()
    @named TY2 = Port()

    #=
    eqs = [
           D(x_1) ~ dx_1
           D(x_2) ~ dx_2
           D(dx_1) ~ ddx_1
           D(dx_2) ~ ddx_2
           D(y_1) ~ dy_1
           D(y_2) ~ dy_2
           D(dy_1) ~ ddy_1
           D(dy_2) ~ ddy_2
           D(A) ~ dA
           D(dA) ~ ddA
           m * (ddx_1 + ddx_2) / 2 ~ fx1 - fx2;
           m * (ddy_1 + ddy_2) / 2 ~ fy1 - fy2 + m * g;
           I * ddA ~ -fy2 * (x_2 - x_1) + fx2 * (y_2 - y_1) + m * g * (x_2 - x_1) / 2;
           x_2 ~ x_1 + l * cos(A);
           y_2 ~ y_1 + l * sin(A);
           dx_1 ~ TX1.v;
           dy_1 ~ TY1.v;
           dx_2 ~ TX2.v;
           dy_2 ~ TY2.v;
           fx1 ~ TX1.f
           fy1 ~ TY1.f
           fx2 ~ -TX2.f
           fy2 ~ -TY2.f]
    =#

    #=
    eqs = [
           m * ((D^2)(x_1) + (D^2)(x_2)) / 2 ~ fx1 - fx2
           m * ((D^2)(y_1) + (D^2)(y_2)) / 2 ~ fy1 - fy2 + m * g
           I * (D^2)(A) ~ -fy2 * (x_2 - x_1) + fx2 * (y_2 - y_1) + m * g * (x_2 - x_1) / 2
           x_2 ~ x_1 + l * cos(A)
           y_2 ~ y_1 + l * sin(A)
           D(x_1) ~ TX1.v
           D(y_1) ~ TY1.v
           D(x_2) ~ TX2.v
           D(y_2) ~ TY2.v
           fx1 ~ TX1.f
           fy1 ~ TY1.f
           fx2 ~ -TX2.f
           fy2 ~ -TY2.f
        ]
    =#

    eqs = [
           m * ((D^2)(x_1) + (D^2)(x_2)) / 2 ~ fx1 - fx2
           m * ((D^2)(y_1) + (D^2)(y_2)) / 2 ~ fy1 - fy2 + m * g
           I * (D^2)(A) ~ -fy2 * (x_2 - x_1) + fx2 * (y_2 - y_1) + m * g * (x_2 - x_1) / 2
           x_2 ~ x_1 + l * cos(A)
           y_2 ~ y_1 + l * sin(A)
           x_1 ~ TX1.s
           y_1 ~ TY1.s
           x_2 ~ TX2.s
           y_2 ~ TY2.s
           fx1 ~ TX1.f
           fy1 ~ TY1.f
           fx2 ~ -TX2.f
           fy2 ~ -TY2.f
        ]

    return ODESystem(eqs, t, vars, pars; name = name, systems = [TX1, TX2, TY1, TY2])
end
