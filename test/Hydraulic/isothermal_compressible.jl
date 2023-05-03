using ModelingToolkit, OrdinaryDiffEq, Test
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC
import ModelingToolkitStandardLibrary.Blocks as B
import ModelingToolkitStandardLibrary.Mechanical.Translational as T

using ModelingToolkitStandardLibrary.Blocks: Parameter

@parameters t
D = Differential(t)

NEWTON = NLNewton(check_div = false, always_new = true, max_iter = 100, relax = 4 // 10)

@testset "Fluid Domain and Tube" begin
    function System(N; bulk_modulus, name)
        pars = @parameters begin bulk_modulus = bulk_modulus end

        systems = @named begin
            fluid = IC.HydraulicFluid(; bulk_modulus)
            stp = B.Step(; height = 10e5, offset = 0, start_time = 0.005, duration = Inf,
                         smooth = true)
            src = IC.Pressure(; p_int = 0)
            vol = IC.FixedVolume(; p_int = 0, vol = 10.0)
        end

        if N == 1
            @named res = IC.TubeBase(; p_int = 0, area = 0.01, length = 500.0)
        else
            @named res = IC.Tube(N; p_int = 0, area = 0.01, length = 500.0)
        end
        push!(systems, res)

        eqs = [connect(stp.output, src.input)
               connect(fluid, src.port)
               connect(src.port, res.port_a)
               connect(res.port_b, vol.port)]

        ODESystem(eqs, t, [], pars; name, systems)
    end

    @named sys1_2 = System(1; bulk_modulus = 2e9)
    @named sys1_1 = System(1; bulk_modulus = 1e9)
    @named sys5_1 = System(5; bulk_modulus = 1e9)

    syss = structural_simplify.([sys1_2, sys1_1, sys5_1])
    probs = [ODEProblem(sys, [], (0, 0.2)) for sys in syss] #ModelingToolkit.missing_variable_defaults(sys)
    sols = [solve(prob, ImplicitEuler(nlsolve = NEWTON); initializealg = NoInit())
            for prob in probs]

    s1_2 = complete(sys1_2)
    s1_1 = complete(sys1_1)
    s5_1 = complete(sys5_1)

    # higher stiffness should compress more quickly and give a higher pressure
    @test sols[1][s1_2.vol.port.p][end] > sols[2][s1_1.vol.port.p][end]

    # N=5 pipe is compressible, will pressurize more slowly
    @test sols[2][s1_1.vol.port.p][end] > sols[3][s5_1.vol.port.p][end]
end

@testset "Valve" begin
    function System(; name)
        pars = []

        systems = @named begin
            fluid = IC.HydraulicFluid()
            sink = IC.FixedPressure(; p = 10e5)
            vol = IC.FixedVolume(; vol = 0.1, p_int = 100e5)
            valve = IC.Valve(; p_a_int = 10e5, p_b_int = 100e5, area_int = 0, Cd = 1e6)
            ramp = B.Ramp(; height = 1, duration = 0.001, offset = 0, start_time = 0.001,
                          smooth = true)
        end

        eqs = [connect(fluid, sink.port)
               connect(sink.port, valve.port_a)
               connect(valve.port_b, vol.port)
               connect(valve.input, ramp.output)]

        ODESystem(eqs, t, [], pars; name, systems)
    end

    @named valve_system = System()
    s = complete(valve_system)
    sys = structural_simplify(valve_system)
    prob = ODEProblem(sys, [], (0, 0.01))
    sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); adaptive = false, dt = 1e-4,
                initializealg = NoInit())

    # the volume should discharge to 10bar
    @test sol[s.vol.port.p][end] ≈ 10e5
end

@testset "DynamicVolume and minimum_volume feature" begin
    function System(; name)
        pars = []

        # DynamicVolume values
        area = 0.01
        length = 0.1

        systems = @named begin
            fluid = IC.HydraulicFluid()
            src1 = IC.Pressure(; p_int = 10e5)
            src2 = IC.Pressure(; p_int = 10e5)

            vol1 = IC.DynamicVolume(+1; p_int = 10e5, area,
                                    dead_volume = 2e-4,
                                    minimum_volume = 2e-4)
            vol2 = IC.DynamicVolume(-1; p_int = 10e5, area,
                                    dead_volume = 2e-4,
                                    minimum_volume = 2e-4, x_int = length)

            mass = T.Mass(; m = 100)

            sin1 = B.Sine(; frequency = 0.5, amplitude = +1e5, offset = 10e5)
            sin2 = B.Sine(; frequency = 0.5, amplitude = -1e5, offset = 10e5)
        end

        eqs = [connect(fluid, src1.port, src2.port)
               connect(src1.port, vol1.port)
               connect(src2.port, vol2.port)
               connect(vol1.flange, mass.flange, vol2.flange)
               connect(src1.input, sin1.output)
               connect(src2.input, sin2.output)]

        ODESystem(eqs, t, [], pars; name, systems)
    end

    @named min_vol_system = System()
    s = complete(min_vol_system)
    sys = structural_simplify(min_vol_system)
    prob = ODEProblem(sys, [], (0, 5.0))

    sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); adaptive = false, dt = 1e-4,
                initializealg = NoInit())

    # begin
    #     fig = Figure()

    #     ax = Axis(fig[1,1], ylabel="position [m]", xlabel="time [s]")
    #     lines!(ax, sol.t, sol[s.vol1.x]; label="vol1")
    #     lines!(ax, sol.t, sol[s.vol2.x]; label="vol2")

    #     ax = Axis(fig[2,1], ylabel="pressure [bar]", xlabel="time [s]")
    #     lines!(ax, sol.t, sol[s.vol1.port.p]/1e5; label="vol1")
    #     lines!(ax, sol.t, sol[s.vol2.port.p]/1e5; label="vol2")
    #     Legend(fig[2,2], ax)
    #     fig    
    # end

    # volume/mass should stop moving at opposite ends
    @test round(sol[s.vol1.vol.x][1]; digits = 2) == 0.0
    @test round(sol[s.vol2.vol.x][1]; digits = 2) == 0.1
    @test round(sol[s.vol1.vol.x][end]; digits = 2) == 0.1
    @test round(sol[s.vol2.vol.x][end]; digits = 2) == 0.0
end

@testset "Actuator System" begin
    function System(use_input, f; name)
        @parameters t

        pars = @parameters begin
            p_s = 285e5
            p_r = 5e5

            A_1 = 908e-4
            A_2 = 360e-4

            p_1 = 17e5
            p_2 = p_1 * A_1 / A_2

            l_1 = 0.7
            l_2 = 2.7
            m_f = 276
            g = 0
            x_f_int = 0

            d = 230e-3

            Cd = 70e5 * 2 * 876 / (876 * 50 / (25e-3 * 2π * 230e-3))^2   #50_000 lpm

            m_piston = 880
            m_body = 1500
        end

        vars = @variables begin ddx(t) = 0 end

        systems = @named begin
            src = IC.FixedPressure(; p = p_s)
            valve = IC.SpoolValve2Way(; p_s_int = p_s, p_a_int = p_1, p_b_int = p_2,
                                      p_r_int = p_r, g, m = m_f, x_int = x_f_int, d, Cd)
            piston = IC.Actuator(; p_a_int = p_1, p_b_int = p_2, area_a = A_1, area_b = A_2,
                                 length_a_int = l_1, length_b_int = l_2, m = m_piston,
                                 g = 0,
                                 x_int = 0, minimum_volume_a = 0, minimum_volume_b = 0)
            body = T.Mass(; m = m_body)
            pipe = IC.Tube(5; p_int = p_2, area = A_2, length = 2.0)
            snk = IC.FixedPressure(; p = p_r)
            pos = T.Position(; s_0 = x_f_int)

            m1 = IC.FlowDivider(; p_int = p_2, n = 3)
            m2 = IC.FlowDivider(; p_int = p_2, n = 3)

            fluid = IC.HydraulicFluid()
        end

        if use_input
            @named input = B.SampledData(Float64)
        else
            @named input = B.TimeVaryingFunction(f; t)
        end

        push!(systems, input)

        eqs = [connect(input.output, pos.input)
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
    prob = ODEProblem(sys, [], (0, 1.0); tofloat = false, jac = true)

    # check the fluid domain
    @test Symbol(defs[s.src.port.ρ]) == Symbol(s.fluid.ρ)
    @test Symbol(defs[s.valve.port_s.ρ]) == Symbol(s.fluid.ρ)
    @test Symbol(defs[s.valve.port_a.ρ]) == Symbol(s.fluid.ρ)
    @test Symbol(defs[s.valve.port_b.ρ]) == Symbol(s.fluid.ρ)
    @test Symbol(defs[s.valve.port_r.ρ]) == Symbol(s.fluid.ρ)
    @test Symbol(defs[s.snk.port.ρ]) == Symbol(s.fluid.ρ)

    dt = 1e-4
    time = 0:dt:0.055
    x = @. 0.9 * (time > 0.015) * (time - 0.015)^2 - 25 * (time > 0.02) * (time - 0.02)^3

    defs[s.input.buffer] = Parameter(x, dt)

    p = Parameter.(ModelingToolkit.varmap_to_vars(defs, parameters(sys); tofloat = false))
    prob = remake(prob; p, tspan = (0, time[end]))
    sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); adaptive = false, dt,
                initializealg = NoInit())

    @test sol[sys.ddx][1] == 0.0
    @test maximum(sol[sys.ddx]) > 500
    @test sol[sys.ddx][end] < 100
end
