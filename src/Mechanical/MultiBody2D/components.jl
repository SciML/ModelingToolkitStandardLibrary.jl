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

    @named M1 = RigidBody2DPort()
    @named M2 = RigidBody2DPort()

    # let ----------------------------------
    Δx = (x2 - x1)/2
    Δy = (y2 - y1)/2

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
           I * ddA ~ Δx*fy2 - Δy*fx2 - Δx*fy1 + Δy*fx1 + M1.T_z + M2.T_z
           

           # geometry
           x2 ~ l * cos(A) + x1
           y2 ~ l * sin(A) + y1
   
           x_cm ~ Δx + x1
           y_cm ~ Δy + y1

           M1.f_x ~ fx1
           M1.f_y ~ fy1
           
           M1.dx ~ dx1
           M1.dy ~ dy1
           M1.dA ~ dA
   
           M2.f_x ~ fx2
           M2.f_y ~ fy2
           
           M2.dx ~ dx2
           M2.dy ~ dy2
           M2.dA ~ dA

           ]

    return ODESystem(eqs, t, vars, pars; name = name, systems = [M1, M2],
    defaults = [M1.dx => 0, M1.dy => 0, M1.dA => 0, M2.dx => 0, M2.dy => 0, M2.dA => 0])
end


function RevoluteJoint(; name, d=0)
    pars = @parameters begin
        d=d
    end

    vars = @variables begin
        dA(t) = 0
        T_z(t) = 0
    end

    @named M1 = RigidBody2DPort()
    @named M2 = RigidBody2DPort()

    eqs = [dA ~ M1.dA - M2.dA
        
        T_z ~ dA*d
        
        M1.T_z ~ +T_z
        M2.T_z ~ -T_z
        
        M1.f_x ~ -M2.f_x
        M1.f_y ~ -M2.f_y
        M1.dx ~ M2.dx
        M1.dy ~ M2.dy]

    return ODESystem(eqs, t, vars, pars; name = name, systems = [M1, M2],
                     defaults = [M1.dx => 0, M1.dy => 0, M1.dA => 0, M2.dx => 0, M2.dy => 0, M2.dA => 0])
end

function MultiBody2Translational(; name)

    @named M = RigidBody2DPort()
    @named T = MechanicalPort()

    eqs = [
            M.f_x ~ -T.f
            
            M.dx ~ T.v
            M.dy ~ 0
            M.dA ~ 0
           ]
    return compose(ODESystem(eqs, t, [], []; name = name, defaults = [M.dx => 0, M.dy => 0, M.dA => 0, T.v => 0, T.f => 0]),
                   M, T)

end

# function FixedFrame(; name)
#     @named T = RigidBody2DPort()
#     eqs = [T.v_x ~ 0
#            T.v_y ~ 0
#            T.dA ~ 0]
#     return compose(ODESystem(eqs, t, [], []; name = name, defaults = [T.v_x => 0, T.v_y => 0, T.dA => 0]),
#                    T)
# end

