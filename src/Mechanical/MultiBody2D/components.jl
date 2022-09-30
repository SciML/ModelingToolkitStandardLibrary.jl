function Link(; name, m, l, I, g)
    
    pars = @parameters begin
        m = m
        l = l
        I = I
        g = g
    end

    vars = @variables begin
        A(t) = 0
        dA(t) = 0
        ddA(t) = 0

        x_cm(t) = l/2
        y_cm(t) = 0

        dx_cm(t) = 0
        dy_cm(t) = 0

        ddx_cm(t) = 0
        ddy_cm(t) = 0

        fx1(t) = 0
        fy1(t) = 0

        x1(t) = 0
        dx1(t) = 0
        
    end

    
    @named TX1 = MechanicalPort()
    

    eqs = [

        D(x_cm) ~ dx_cm
        D(y_cm) ~ dy_cm

        D(dx_cm) ~ ddx_cm
        D(dy_cm) ~ ddy_cm

        D(A) ~ dA
        D(dA) ~ ddA

        # x forces
        m*ddx_cm ~ fx1

        # y forces
        m*ddy_cm ~ m*g + fy1

        # torques
        I * ddA ~ m * g * (l/2) * cos(A)

        # geometry
        x_cm ~ l * cos(A)/2 + x1
        y_cm ~ l * sin(A)/2
        
        TX1.f ~ fx1
        D(x1) ~ TX1.v
        
    ]


    return ODESystem(eqs, t, vars, pars; name = name, systems = [TX1], defaults = [TX1.v => 0])
end
