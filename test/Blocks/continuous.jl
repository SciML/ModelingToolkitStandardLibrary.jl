using ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks

@parameters t

#=
Testing strategy:
The general strategy is to test systems using simple intputs where the solution
is known on closed form. For algebraic systems (without differential variables),
an integrator with a constant input is often used together with the system under test. 
=#

@testset "Constant" begin
    @named c = Constant(; k=1)
    @named int = Integrator()
    @named iosys = ODESystem(connect(c.output, int.input), t, systems=[int, c])
    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[int.x=>1.0], (0.0, 1.0))
    sol = solve(prob, Rodas4())
    @test all(sol[int.output.u] .≈ 1)
end

@testset "Derivative" begin
    @named source = Sine(; frequency=1)
    @named int = Integrator(; k=1)
    @named der = Derivative(; k=1, T=0.001)
    @named iosys = ODESystem([
        connect(source.output, der.input),
        connect(der.output, int.input),
        ],
        t,
        systems=[int, source, der],
    )
    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())
    @test all(isapprox.(sol[source.output.u], sol[int.output.u], atol=1e-1))
end

@testset "PT1" begin
    pt1_func(t, k, T) = k * (1 - exp(-t / T)) # Known solution to first-order system
    
    k, T = 1.0, 0.1
    @named c = Constant(; k=1)
    @named pt1 = FirstOrder(; k=k, T=T)
    @named iosys = ODESystem(connect(c.output, pt1.input), t, systems=[pt1, c])
    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test sol[pt1.output.u] ≈ pt1_func.(sol.t, k, T) atol=1e-3
end

@testset "PT2" begin
    # Known solution to second-order system
    function pt2_func(t, k, w, d)
        y = if d==0
            -k*(-1 + cos(t*w))
        else
            d = complex(d)
            real(k*(1 + (-cosh(sqrt(-1 + d^2)*t*w) - (d*sinh(sqrt(-1 + d^2)*t*w))/sqrt(-1 + d^2))/exp(d*t*w)))
        end
    end

    k, w, d = 1.0, 1.0, 0.5
    @named c = Constant(; k=1)
    @named pt2 = SecondOrder(; k=k, w=w, d=d)
    @named iosys = ODESystem(connect(c.output, pt2.input), t, systems=[pt2, c])
    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test sol[pt2.output.u] ≈ pt2_func.(sol.t, k, w, d) atol=1e-3
end

@testset "StateSpace" begin
    A = [0 1;-1 -0.5]
    B = [0, 1]
    C = [0.9 1;]
    D = [0;;]
    @named ss = StateSpace(;A,B,C,D,x_start=zeros(2))
    @named c = Constant(; k=1)
    @named model = ODESystem(
        [
           connect(c.output, ss.input),
        ], 
        t, 
        systems=[ss, c]
    )
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    # initial condition
    @test sol[ss.x[1]][1] ≈ 0 atol=1e-3
    @test sol[ss.x[2]][1] ≈ 0 atol=1e-3
    # equilibrium point is at [1, 0]
    @test sol[ss.x[1]][end] ≈ 1 atol=1e-3
    @test sol[ss.x[2]][end] ≈ 0 atol=1e-3
end

"""Second order demo plant"""
function Plant(;name, x_start=zeros(2))
    @named input = RealInput()
    @named output = RealOutput()
    D = Differential(t)
    sts = @variables x1(t)=x_start[1] x2(t)=x_start[2]
    eqs= [
        D(x1) ~ x2
        D(x2) ~ -x1 - 0.5 * x2 + input.u
        output.u ~ 0.9 * x1 + x2
    ]
    compose(ODESystem(eqs, t, sts, []; name), [input, output])
end

@testset "PI" begin
    re_val = 2
    @named ref = Constant(; k=re_val)
    @named pi_controller = PI(k=1, T=1)
    @named plant = Plant()
    @named fb = Feedback()
    @named model = ODESystem(
        [
            connect(ref.output, fb.input1), 
            connect(plant.output, fb.input2),
            connect(fb.output, pi_controller.err_input), 
            connect(pi_controller.ctr_output, plant.input), 
        ], 
        t, 
        systems=[pi_controller, plant, ref, fb]
    )
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
    @test sol[plant.output.u][end] ≈ re_val atol=1e-3 # zero control error after 100s
end

@testset "PID" begin
    re_val = 2
    @named ref = Constant(; k=re_val)
    @named pid_controller = PID(k=3, Ti=0.5, Td=100)
    @named plant = Plant()
    @named fb = Feedback()
    @named model = ODESystem(
        [
            connect(ref.output, fb.input1), 
            connect(plant.output, fb.input2),
            connect(fb.output, pid_controller.err_input), 
            connect(pid_controller.ctr_output, plant.input), 
        ], 
        t, 
        systems=[pid_controller, plant, ref, fb]
    )
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())
    @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
    @test sol[plant.output.u][end] ≈ re_val atol=1e-3 # zero control error after 100s

    @testset "PI" begin
        @named pid_controller = PID(k=3, Ti=0.5, Td=false)
        @named model = ODESystem(
            [
                connect(ref.output, fb.input1), 
                connect(plant.output, fb.input2),
                connect(fb.output, pid_controller.err_input), 
                connect(pid_controller.ctr_output, plant.input), 
            ], 
            t, 
            systems=[pid_controller, plant, ref, fb]
        )
        sys = structural_simplify(model)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())
        @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
        @test sol[plant.output.u][end] ≈ re_val atol=1e-3 # zero control error after 100s
    end

    @testset "PD" begin
        @named pid_controller = PID(k=10, Ti=false, Td=1)
        @named model = ODESystem(
            [
                connect(ref.output, fb.input1), 
                connect(plant.output, fb.input2),
                connect(fb.output, pid_controller.err_input), 
                connect(pid_controller.ctr_output, plant.input), 
            ], 
            t, 
            systems=[pid_controller, plant, ref, fb]
        )
        sys = structural_simplify(model)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())
        @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
        @test sol[plant.output.u][end] > 1 # without I there will be a steady-state error
    end
end

@test_skip begin 
    k = 2
    Ti = 0.5
    Td = 0.7
    wp = 1
    wd = 1
    Ni = √(Td / Ti)
    Nd = 12
    y_max = Inf
    y_min = -Inf
    u_r = sin(t)
    u_y = 0
    function solve_with_input(; u_r, u_y, 
        controller = PID(; k, Ti, Td, wp, wd, Ni, Nd, y_max, y_min, name=:controller)
    )
        @test count(ModelingToolkit.isinput, states(controller)) == 5 # 2 in PID, 1 sat, 1 I, 1 D
        @test count(ModelingToolkit.isoutput, states(controller)) == 4
        # TODO: check number of unbound inputs when available, should be 2
        @named iosys = ODESystem([controller.u_r~u_r, controller.u_y~u_y], t, systems=[controller])
        sys = structural_simplify(iosys)
        prob = ODEProblem(sys, Pair[], (0.0, 10.0))
        sol = solve(prob, Rodas4(), saveat=0:0.1:10)
        controller, sys, sol
    end

    # linearity in u_r
    controller, sys, sol1 = solve_with_input(u_r=sin(t), u_y=0)
    controller, sys, sol2 = solve_with_input(u_r=2sin(t), u_y=0)
    @test sum(abs, sol1[controller.ea]) < eps() # This is the acutator model error due to saturation
    @test 2sol1[controller.y] ≈ sol2[controller.y] rtol=1e-3 # linearity in u_r

    # linearity in u_y
    controller, sys, sol1 = solve_with_input(u_y=sin(t), u_r=0)
    controller, sys, sol2 = solve_with_input(u_y=2sin(t), u_r=0)
    @test sum(abs, sol1[controller.ea]) < eps() # This is the acutator model error due to saturation
    @test 2sol1[controller.y] ≈ sol2[controller.y] rtol=1e-3 # linearity in u_y

    # zero error
    controller, sys, sol1 = solve_with_input(u_y=sin(t), u_r=sin(t))
    @test sum(abs, sol1[controller.y]) ≈ 0 atol=sqrt(eps()) 

    # test saturation
    controller, sys, sol1 = solve_with_input(; u_r=10sin(t), u_y=0, 
        controller = PID(; k, Ti, Td, wp, wd=0, Ni, Nd, y_max=10, y_min=-10, name=:controller)
    )
    @test extrema(sol1[controller.y]) == (-10, 10)


    # test P set-point weighting
    controller, sys, sol1 = solve_with_input(; u_r=sin(t), u_y=0, 
        controller = PID(; k, Ti, Td, wp=0, wd, Ni, Nd, y_max, y_min, name=:controller)
    )
    @test sum(abs, sol1[controller.ep]) ≈ 0 atol=sqrt(eps()) 

    # test D set-point weighting
    controller, sys, sol1 = solve_with_input(; u_r=sin(t), u_y=0, 
        controller = PID(; k, Ti, Td, wp, wd=0, Ni, Nd, y_max, y_min, name=:controller)
    )
    @test sum(abs, sol1[controller.ed]) ≈ 0 atol=sqrt(eps()) 


    # zero integral gain
    controller, sys, sol1 = solve_with_input(; u_r=sin(t), u_y=0, 
        controller = PID(; k, Ti=false, Td, wp, wd, Ni, Nd, y_max, y_min, name=:controller)
    )
    @test isapprox(sum(abs, sol1[controller.I.y]), 0, atol=sqrt(eps()))
    

    # zero derivative gain
    @test_skip begin # During the resolution of the non-linear system, the evaluation of the following equation(s) resulted in a non-finite number: [5]
        controller, sys, sol1 = solve_with_input(; u_r=sin(t), u_y=0, 
            controller = PID(; k, Ti, Td=false, wp, wd, Ni, Nd, y_max, y_min, name=:controller)
        )
        @test isapprox(sum(abs, sol1[controller.D.y]), 0, atol=sqrt(eps()))
    end

    # Tests below can be activated when the concept of unbound_inputs exists in MTK
    # @test isequal(Set(unbound_inputs(controller)), @nonamespace(Set([controller.u_r, controller.u_y])))
    # @test isempty(unbound_inputs(sys))
    # @test isequal(bound_inputs(sys), inputs(sys))
    # @test isequal(
    #     Set(bound_inputs(sys)),
    #     Set([controller.u_r, controller.u_y, controller.I.u, controller.D.u, controller.sat.u])
    #     )
end

@testset "LimPI" begin
    re_val = 1
    @named ref = Constant(; k=re_val)
    @named pi_controller_lim = LimPI(k=3, T=0.5, u_max=1.5, u_min=-1.5, Ta=0.1)
    @named pi_controller = PI(k=3, T=0.5)
    @named sat = Limiter(y_max=1.5, y_min=-1.5)
    @named plant = Plant()
    @named fb = Feedback()

    # without anti-windup measure
    sol = let
        @named model = ODESystem(
            [
                connect(ref.output, fb.input1), 
                connect(plant.output, fb.input2),
                connect(fb.output, pi_controller.err_input), 
                connect(pi_controller.ctr_output, sat.input), 
                connect(sat.output, plant.input), 
            ], 
            t, 
            systems=[pi_controller, plant, ref, fb, sat]
        )
        sys = structural_simplify(model)
        prob = ODEProblem(sys, Pair[], (0.0, 20.0))
        sol = solve(prob, Rodas4())
    end

    # with anti-windup measure
    sol_lim = let
        @named model = ODESystem(
            [
                connect(ref.output, fb.input1), 
                connect(plant.output, fb.input2),
                connect(fb.output, pi_controller_lim.err_input), 
                connect(pi_controller_lim.ctr_output, sat.input), 
                connect(sat.output, plant.input), 
            ], 
            t, 
            systems=[pi_controller_lim, plant, ref, fb, sat]
        )
        sys = structural_simplify(model)
        prob = ODEProblem(sys, Pair[], (0.0, 20.0))
        sol = solve(prob, Rodas4())
    end

    @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
    @test all(isapprox.(sol_lim[ref.output.u], re_val, atol=1e-3))  # check reference
    @test sol[plant.output.u][end] ≈ re_val atol=1e-3 # zero control error after 100s
    @test sol_lim[plant.output.u][end] ≈ re_val atol=1e-3 # zero control error after 100s
    @test all(-1.5 .<= sol_lim[pi_controller_lim.ctr_output.u] .<= 1.5) # test limit

    # Plots.plot(sol; vars=[plant.output.u]) # without anti-windup measure
    # Plots.plot!(sol_lim; vars=[plant.output.u]) # with anti-windup measure
end

@testset "LimPID" begin
    re_val = 1
    @named ref = Constant(; k=re_val)
    @named pid_controller = LimPID(k=3, Ti=0.5, Td=100, u_max=1.5, u_min=-1.5, Ni=0.1/0.5)
    @named plant = Plant()
    @named model = ODESystem(
        [
            connect(ref.output, pid_controller.reference), 
            connect(plant.output, pid_controller.measurement),
            connect(pid_controller.ctr_output, plant.input), 
        ], 
        t, 
        systems=[pid_controller, plant, ref]
    )
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0.0, 100.0))
    sol = solve(prob, Rodas4())

    # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
    @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
    @test sol[plant.output.u][end] ≈ re_val atol=1e-3 # zero control error after 100s
    @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit

    @testset "PI" begin
        @named pid_controller = LimPID(k=3, Ti=0.5, Td=false, u_max=1.5, u_min=-1.5, Ni=0.1/0.5)
        @named model = ODESystem(
            [
                connect(ref.output, pid_controller.reference), 
                connect(plant.output, pid_controller.measurement),
                connect(pid_controller.ctr_output, plant.input), 
            ], 
            t, 
            systems=[pid_controller, plant, ref]
        )
        sys = structural_simplify(model)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())

        # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
        @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
        @test sol[plant.output.u][end] ≈ re_val atol=1e-3 # zero control error after 100s
        @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit
    end
    @testset "PD" begin
        @named pid_controller = LimPID(k=10, Ti=false, Td=1, u_max=1.5, u_min=-1.5)
        @named model = ODESystem(
            [
                connect(ref.output, pid_controller.reference), 
                connect(plant.output, pid_controller.measurement),
                connect(pid_controller.ctr_output, plant.input), 
            ], 
            t, 
            systems=[pid_controller, plant, ref]
        )
        sys = structural_simplify(model)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())

        # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
        @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
        @test sol[plant.output.u][end] > 0.5 # without I there will be a steady-state error
        @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit
    end
    @testset "set-point weights" begin
        @testset "wp" begin
            @named pid_controller = LimPID(k=3, Ti=0.5, Td=100, u_max=1.5, u_min=-1.5, Ni=0.1/0.5, wp=0, wd=1)
            @named model = ODESystem(
                [
                    connect(ref.output, pid_controller.reference), 
                    connect(plant.output, pid_controller.measurement),
                    connect(pid_controller.ctr_output, plant.input), 
                ], 
                t, 
                systems=[pid_controller, plant, ref]
            )
            sys = structural_simplify(model)
            prob = ODEProblem(sys, Pair[], (0.0, 100.0))
            sol = solve(prob, Rodas4())

            # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
            @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
            sol[pid_controller.addP.output.u] == -sol[pid_controller.measurement.u]
            @test sol[plant.output.u][end] ≈ re_val atol=1e-3 # zero control error after 100s
            @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit
        end
        @testset "wd" begin
            @named pid_controller = LimPID(k=3, Ti=0.5, Td=100, u_max=1.5, u_min=-1.5, Ni=0.1/0.5, wp=1, wd=0)
            @named model = ODESystem(
                [
                    connect(ref.output, pid_controller.reference), 
                    connect(plant.output, pid_controller.measurement),
                    connect(pid_controller.ctr_output, plant.input), 
                ], 
                t, 
                systems=[pid_controller, plant, ref]
            )
            sys = structural_simplify(model)
            prob = ODEProblem(sys, Pair[], (0.0, 100.0))
            sol = solve(prob, Rodas4())

            # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
            @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
            @test sol[plant.output.u][end] ≈ re_val atol=1e-3 # zero control error after 100s
            sol[pid_controller.addD.output.u] == -sol[pid_controller.measurement.u]
            @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit
        end
    end
    @testset "PI without AWM" begin
        @named pid_controller = LimPID(k=3, Ti=0.5, Td=false, u_max=1.5, u_min=-1.5, Ni=Inf)
        @named model = ODESystem(
            [
                connect(ref.output, pid_controller.reference), 
                connect(plant.output, pid_controller.measurement),
                connect(pid_controller.ctr_output, plant.input), 
            ], 
            t, 
            systems=[pid_controller, plant, ref]
        )
        sys = structural_simplify(model)
        prob = ODEProblem(sys, Pair[], (0.0, 100.0))
        sol = solve(prob, Rodas4())

        # Plots.plot(sol, vars=[plant.output.u, plant.input.u])
        @test all(isapprox.(sol[ref.output.u], re_val, atol=1e-3))  # check reference
        @test sol[plant.output.u][end] ≈ re_val atol=1e-3 # zero control error after 100s
        @test all(-1.5 .<= sol[pid_controller.ctr_output.u] .<= 1.5) # test limit
    end
end