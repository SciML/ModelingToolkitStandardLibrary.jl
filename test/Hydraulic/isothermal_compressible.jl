using ModelingToolkit, OrdinaryDiffEq, Test
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B
import ModelingToolkitStandardLibrary.Mechanical.Translational as T

using ModelingToolkitStandardLibrary.Blocks: Parameter

@parameters t
D = Differential(t)

NEWTON = NLNewton(check_div = false, always_new = true, max_iter = 100, relax = 9 // 10)

@testset "Fluid Domain and Tube" begin
    function System(N; bulk_modulus, name)
        pars = @parameters begin
            bulk_modulus = bulk_modulus
        end

        systems = @named begin
            fluid = IC.HydraulicFluid(; bulk_modulus)
            stp = B.Step(; height = 10e5, offset = 0, start_time = 0.005, duration = Inf,
                smooth = true)
            src = IC.Pressure(; p_int = 0)
            vol = IC.FixedVolume(; p_int = 0, vol = 10.0)
            res = IC.Tube(N; p_int = 0, area = 0.01, length = 50.0)
        end

        eqs = [connect(stp.output, src.p)
            connect(fluid, src.port)
            connect(src.port, res.port_a)
            connect(res.port_b, vol.port)]

        ODESystem(eqs, t, [], pars; name, systems)
    end

    @named sys1_2 = System(1; bulk_modulus = 2e9)
    @named sys1_1 = System(1; bulk_modulus = 1e9)
    @named sys5_1 = System(5; bulk_modulus = 1e9)

    syss = structural_simplify.([sys1_2, sys1_1, sys5_1])
    probs = [ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 0.05))
             for sys in syss] #
    sols = [solve(prob, ImplicitEuler(nlsolve = NEWTON); initializealg = NoInit(),
        dt = 1e-4, adaptive = false)
            for prob in probs]

    s1_2 = complete(sys1_2)
    s1_1 = complete(sys1_1)
    s5_1 = complete(sys5_1)

    # higher stiffness should compress more quickly and give a higher pressure
    @test sols[1][s1_2.vol.port.p][end] > sols[2][s1_1.vol.port.p][end]

    # N=5 pipe is compressible, will pressurize more slowly
    @test sols[2][s1_1.vol.port.p][end] > sols[3][s5_1.vol.port.p][end]

    #=
    fig = Figure()
    ax = Axis(fig[1,1])
    # hlines!(ax, 10e5)
    lines!(ax, sols[1][s1_2.vol.port.p])
    lines!(ax, sols[2][s1_1.vol.port.p])
    lines!(ax, sols[3][s5_1.vol.port.p])
    fig
    =#

end

@testset "Valve" begin
    function System(; name)
        pars = []

        systems = @named begin
            fluid = IC.HydraulicFluid()
            sink = IC.FixedPressure(; p = 10e5)
            vol = IC.FixedVolume(; vol = 0.1, p_int = 100e5)
            valve = IC.Valve(; p_a_int = 10e5, p_b_int = 100e5, area_int = 0, Cd = 1e5,
                minimum_area = 0)
            ramp = B.Ramp(; height = 1, duration = 0.001, offset = 0, start_time = 0.001,
                smooth = true)
        end

        eqs = [connect(fluid, sink.port)
            connect(sink.port, valve.port_a)
            connect(valve.port_b, vol.port)
            connect(valve.area, ramp.output)]

        ODESystem(eqs, t, [], pars; name, systems)
    end

    @named valve_system = System()
    s = complete(valve_system)
    sys = structural_simplify(valve_system)
    prob = ODEProblem(sys, [], (0, 0.01))
    sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); adaptive = false, dt = 1e-4,
        initializealg = NoInit())

    # the volume should discharge to 10bar
    @test sol[s.vol.port.p][end]≈10e5 atol=1e5
end

@testset "Volume" begin
    dt = 1e-4 #s
    t_end = 0.2 #s
    time = 0:dt:t_end

    using DataInterpolations
    u = include("ref/dm.jl")
    dm_fun = LinearInterpolation(u, time)

    dx = 4.71238898038469
    drho = -30556.685668435468
    dm = dm_fun(0.0)

    function MassVolume(; name, dx, drho, dm)
        pars = @parameters begin
            A = 0.01 #m²
            x₀ = 1.0 #m
            M = 10_000 #kg
            g = 9.807 #m/s²
            amp = 5e-2 #m
            f = 15 #Hz   
            p_int = M * g / A
            dx = dx
            drho = drho
            dm = dm
        end
        vars = []
        systems = @named begin
            fluid = IC.HydraulicFluid(; density = 876, bulk_modulus = 1.2e9)
            mass = T.Mass(; v = dx, m = M, g = -g)
            vol = IC.Volume(; area = A, x = x₀, p = p_int, dx, drho, dm)
            mass_flow = IC.MassFlow(; p_int)
            mass_flow_input = B.TimeVaryingFunction(; f = dm_fun)
        end

        eqs = [connect(mass.flange, vol.flange)
            connect(vol.port, mass_flow.port)
            connect(mass_flow.dm, mass_flow_input.output)
            connect(mass_flow.port, fluid)]

        return ODESystem(eqs, t, vars, pars; systems, name)
    end

    @named odesys = MassVolume(; dx, drho, dm)
    sys = structural_simplify(odesys) |> complete
    prob = ODEProblem(sys, [], (0, 0.2))
    sol = solve(prob, Rodas5P(); dt = 1e-4, adaptive = false)

    correct_x = 5e-2 * sin.(2π * sol.t * 15) .+ 1.0
    err = sol[sys.vol.x] .- correct_x
    @test OrdinaryDiffEq.norm(err) < 1e-4

    #=
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, sol.t, sol[sys.vol.x]; linewidth=2)
    lines!(ax, sol.t, correct_x)
    fig
    =#

end

@testset "DynamicVolume and minimum_volume feature" begin
    function System(N; name, area = 0.01, length = 0.1, damping_volume)
        pars = []

        # DynamicVolume values
        systems = @named begin
            fluid = IC.HydraulicFluid(; bulk_modulus = 1e9)
            src1 = IC.Pressure(; p_int = 10e5)
            src2 = IC.Pressure(; p_int = 10e5)

            vol1 = IC.DynamicVolume(N; direction = +1,
                p_int = 10e5, area,
                x_int = length,
                x_max = length * 2,
                x_min = length * 0.1,
                x_damp = damping_volume / area + length * 0.1)
            vol2 = IC.DynamicVolume(N; direction = -1,
                p_int = 10e5, area,
                x_int = length,
                x_max = length * 2,
                x_min = length * 0.1,
                x_damp = damping_volume / area + length * 0.1)

            mass = T.Mass(; m = 10)

            sin1 = B.Sine(; frequency = 0.5, amplitude = +0.5e5, offset = 10e5)
            sin2 = B.Sine(; frequency = 0.5, amplitude = -0.5e5, offset = 10e5)
        end

        eqs = [connect(fluid, src1.port)
            connect(fluid, src2.port)
            connect(src1.port, vol1.port)
            connect(src2.port, vol2.port)
            connect(vol1.flange, mass.flange, vol2.flange)
            connect(src1.p, sin1.output)
            connect(src2.p, sin2.output)]

        ODESystem(eqs, t, [], pars; name, systems)
    end

    for N in [1, 2]
        for damping_volume in [0.01 * 0.1 * 0.25]
            @named system = System(N; damping_volume)
            s = complete(system)
            sys = structural_simplify(system)
            prob = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 5),
                [s.vol1.Cd_reverse => 0.1, s.vol2.Cd_reverse => 0.1];
                jac = true)

            @time sol = solve(prob,
                ImplicitEuler(nlsolve = NLNewton(check_div = false,
                    always_new = true,
                    max_iter = 10,
                    relax = 9 // 10));
                dt = 0.0001, adaptive = false, initializealg = NoInit())

            # begin
            #     fig = Figure()

            #     ax = Axis(fig[1,1], ylabel="position [m]", xlabel="time [s]")
            #     lines!(ax, sol.t, sol[s.vol1.x]; label="vol1")
            #     lines!(ax, sol.t, sol[s.vol2.x]; label="vol2")
            #     Legend(fig[1,2], ax)

            #     ax = Axis(fig[2,1], ylabel="pressure [bar]", xlabel="time [s]")
            #     lines!(ax, sol.t, sol[s.vol1.damper.port_a.p]/1e5; label="vol1")
            #     lines!(ax, sol.t, sol[s.vol2.damper.port_a.p]/1e5; label="vol2")
            #     ylims!(ax, 10-0.5, 10+0.5)

            #     ax = Axis(fig[3,1], ylabel="area", xlabel="time [s]")
            #     lines!(ax, sol.t, sol[s.vol1.damper.area]; label="area 1")
            #     lines!(ax, sol.t, sol[s.vol2.damper.area]; label="area 2")

            #     display(fig)
            # end

            i1 = round(Int, 1 / 1e-4)
            i2 = round(Int, 2 / 1e-4)
            i3 = round(Int, 3 / 1e-4)
            i4 = round(Int, 4 / 1e-4)

            # volume/mass should stop moving at opposite ends
            @test round(sol[s.vol1.x][1]; digits = 2) == 0.1
            @test round(sol[s.vol2.x][1]; digits = 2) == 0.1
            @test round(sol[s.vol1.x][i1]; digits = 2) == +0.19
            @test round(sol[s.vol2.x][i1]; digits = 2) == +0.01
            @test round(sol[s.vol1.x][i2]; digits = 2) == +0.01
            @test round(sol[s.vol2.x][i2]; digits = 2) == +0.19
            @test round(sol[s.vol1.x][i3]; digits = 2) == +0.19
            @test round(sol[s.vol2.x][i3]; digits = 2) == +0.01
            @test round(sol[s.vol1.x][i4]; digits = 2) == +0.01
            @test round(sol[s.vol2.x][i4]; digits = 2) == +0.19
        end
    end
end

@testset "Actuator System" begin
    function System(use_input, f; name)
        @parameters t

        pars = @parameters begin
            p_s = 200e5
            p_r = 5e5

            A_1 = 360e-4
            A_2 = 360e-4

            p_1 = 45e5
            p_2 = 45e5

            l_1 = 0.01
            l_2 = 0.05
            m_f = 250
            g = 0

            d = 100e-3

            Cd = 0.01

            m_piston = 880
        end

        vars = @variables begin
            ddx(t) = 0
        end

        systems = @named begin
            src = IC.FixedPressure(; p = p_s)
            valve = IC.SpoolValve2Way(; p_s_int = p_s, p_a_int = p_1, p_b_int = p_2,
                p_r_int = p_r, g, m = m_f, x_int = 0, d, Cd)
            piston = IC.Actuator(5;
                p_a_int = p_1,
                p_b_int = p_2,
                area_a = A_1,
                area_b = A_2,
                length_a_int = l_1,
                length_b_int = l_2,
                m = m_piston,
                g = 0,
                x_int = 0,
                minimum_volume_a = A_1 * 1e-3,
                minimum_volume_b = A_2 * 1e-3,
                damping_volume_a = A_1 * 5e-3,
                damping_volume_b = A_2 * 5e-3)
            body = T.Mass(; m = 1500)
            pipe = IC.Tube(5; p_int = p_2, area = A_2, length = 2.0)
            snk = IC.FixedPressure(; p = p_r)
            pos = T.Position()

            m1 = IC.FlowDivider(; p_int = p_2, n = 3)
            m2 = IC.FlowDivider(; p_int = p_2, n = 3)

            fluid = IC.HydraulicFluid()
        end

        if use_input
            @named input = B.SampledData(Float64)
        else
            @named input = B.TimeVaryingFunction(f)
        end

        push!(systems, input)

        eqs = [connect(input.output, pos.s)
            connect(valve.flange, pos.flange)
            connect(valve.port_a, piston.port_a)
            connect(piston.flange, body.flange)
            connect(piston.port_b, m1.port_a)
            connect(m1.port_b, pipe.port_b)
            connect(pipe.port_a, m2.port_b)
            connect(m2.port_a, valve.port_b)
            connect(src.port, valve.port_s)
            connect(snk.port, valve.port_r)
            connect(fluid, src.port, snk.port)
            D(body.v) ~ ddx]

        ODESystem(eqs, t, vars, pars; name, systems)
    end

    @named system = System(true, nothing)

    sys = structural_simplify(system)
    defs = ModelingToolkit.defaults(sys)
    s = complete(system)
    dt = 1e-4
    time = 0:dt:0.1
    x = @. 0.9 * (time > 0.015) * (time - 0.015)^2 - 25 * (time > 0.02) * (time - 0.02)^3
    defs[s.input.buffer] = Parameter(x, dt)
    #prob = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 0.1);
    #                  tofloat = false)#, jac = true)
    prob = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 0.1);
        tofloat = false, jac = true)

    # check the fluid domain
    @test Symbol(defs[s.src.port.ρ]) == Symbol(s.fluid.ρ)
    @test Symbol(defs[s.valve.port_s.ρ]) == Symbol(s.fluid.ρ)
    @test Symbol(defs[s.valve.port_a.ρ]) == Symbol(s.fluid.ρ)
    @test Symbol(defs[s.valve.port_b.ρ]) == Symbol(s.fluid.ρ)
    @test Symbol(defs[s.valve.port_r.ρ]) == Symbol(s.fluid.ρ)
    @test Symbol(defs[s.snk.port.ρ]) == Symbol(s.fluid.ρ)

    # defs[s.piston.Cd_reverse] = 0.1

    prob = remake(prob; tspan = (0, time[end]))
    @time sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); adaptive = false, dt,
        initializealg = NoInit())

    @test sol[sys.ddx][1] == 0.0
    @test maximum(sol[sys.ddx]) > 200
    @test sol[s.piston.x][end]≈0.05 atol=0.01
end

@testset "Prevent Negative Pressure" begin
    @component function System(; name)
        pars = @parameters let_gas = 1

        # TODO: DynamicVolume doesn't work with N>1
        systems = @named begin
            fluid = IC.HydraulicFluid(; let_gas)
            vol = IC.DynamicVolume(1; p_int = 100e5, area = 0.001, x_int = 0.05,
                x_max = 0.1, x_damp = 0.02, x_min = 0.01, direction = +1)
            # vol = IC.Volume(; x = 0.05, p=100e5, area = 0.001)
            mass = T.Mass(; m = 100, g = -9.807, s = 0.05)
            cap = IC.Cap(; p_int = 100e5)
        end

        eqs = [connect(fluid, cap.port)
            connect(cap.port, vol.port)
            connect(vol.flange, mass.flange)]

        ODESystem(eqs, t, [], pars; name, systems)
    end

    @named system = System()
    s = complete(system)
    sys = structural_simplify(system)
    prob1 = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 0.05))
    prob2 = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 0.05),
        [s.let_gas => 0])

    @time sol1 = solve(prob1, ImplicitEuler(nlsolve = NEWTON); adaptive = false, dt = 1e-4)
    @time sol2 = solve(prob2, Rodas4())

    # case 1: no negative pressure will only have gravity pulling mass back down
    # case 2: with negative pressure, added force pulling mass back down
    # - case 1 should push the mass higher
    @test maximum(sol1[s.mass.s]) > maximum(sol2[s.mass.s])

    # case 1 should prevent negative pressure less than -1000
    @test minimum(sol1[s.vol.port.p]) > -1000
    @test minimum(sol2[s.vol.port.p]) < -1000

    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, sol1.t, sol1[s.vol.port.p])
    lines!(ax, sol2.t, sol2[s.vol.port.p])

    ax = Axis(fig[1, 2])
    lines!(ax, sol1.t, sol1[s.mass.s])
    lines!(ax, sol2.t, sol2[s.mass.s])
    fig
end

#TODO: Test Valve Inversion
