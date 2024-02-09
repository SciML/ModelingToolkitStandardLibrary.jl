using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Thermal
import ModelingToolkitStandardLibrary
using ModelingToolkit, OrdinaryDiffEq, Test
using OrdinaryDiffEq: ReturnCode.Success
using DynamicQuantities: @u_str
# using Plots

@parameters t [unit = u"s"]
D = Differential(t)

@testset "DC motor" begin
    R = 0.5
    L = 4.5e-3
    k = 0.5
    J = 0.02
    f = 0.01
    V_step = 10
    tau_L_step = -3
    @named ground = Ground()
    @named source = Voltage()
    @named voltage_step = Blocks.Step(height = V_step, start_time = 0, output__unit = u"V")
    @named R1 = Resistor(R = R)
    @named L1 = Inductor(L = L, i = 0.0)
    @named emf = EMF(k = k)
    @named fixed = Fixed()
    @named load = Torque(use_support = false)
    @named load_step = Blocks.Step(height = tau_L_step, start_time = 3, output__unit = u"N*m")
    @named inertia = Inertia(J = J)
    @named friction = Damper(d = f)

    connections = [connect(fixed.flange, emf.support, friction.flange_b)
        connect(emf.flange, friction.flange_a, inertia.flange_a)
        connect(inertia.flange_b, load.flange)
        connect(load_step.output, load.tau)
        connect(voltage_step.output, source.V)
        connect(source.p, R1.p)
        connect(R1.n, L1.p)
        connect(L1.n, emf.p)
        connect(emf.n, source.n, ground.g)]

    @named model = ODESystem(connections, t,
        systems = [
            ground,
            voltage_step,
            source,
            R1,
            L1,
            emf,
            fixed,
            load,
            load_step,
            inertia,
            friction,
        ])
    sys = structural_simplify(model)

    prob = ODEProblem(sys, [], (0, 6.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    # EMF equations
    @test -0.5 .* sol[emf.i] == sol[emf.flange.tau]
    @test sol[emf.v] == 0.5 .* sol[emf.w]
    # test steady-state values
    dc_gain = [f/(k^2 + f * R) k/(k^2 + f * R); k/(k^2 + f * R) -R/(k^2 + f * R)]
    idx_t = findfirst(sol.t .> 2.5)
    @test sol[inertia.w][idx_t]≈(dc_gain * [V_step; 0])[2] rtol=1e-3
    @test sol[emf.i][idx_t]≈(dc_gain * [V_step; 0])[1] rtol=1e-3
    idx_t = findfirst(sol.t .> 5.5)
    @test sol[inertia.w][idx_t]≈(dc_gain * [V_step; -tau_L_step])[2] rtol=1e-3
    @test sol[emf.i][idx_t]≈(dc_gain * [V_step; -tau_L_step])[1] rtol=1e-3

    prob = DAEProblem(sys, D.(unknowns(sys)) .=> 0.0, Pair[], (0, 6.0))
    sol = solve(prob, DFBDF())
    @test sol.retcode == Success
    # EMF equations
    @test -0.5 .* sol[emf.i] == sol[emf.flange.tau]
    @test sol[emf.v] == 0.5 .* sol[emf.w]
    # test steady-state values
    dc_gain = [f/(k^2 + f * R) k/(k^2 + f * R); k/(k^2 + f * R) -R/(k^2 + f * R)]
    idx_t = findfirst(sol.t .> 2.5)
    @test sol[inertia.w][idx_t]≈(dc_gain * [V_step; 0])[2] rtol=1e-3
    @test sol[emf.i][idx_t]≈(dc_gain * [V_step; 0])[1] rtol=1e-3
    idx_t = findfirst(sol.t .> 5.5)
    @test sol[inertia.w][idx_t]≈(dc_gain * [V_step; -tau_L_step])[2] rtol=1e-3
    @test sol[emf.i][idx_t]≈(dc_gain * [V_step; -tau_L_step])[1] rtol=1e-3

    # p1 = Plots.plot(sol, vars=[inertia.w], ylabel="Angular Vel. in rad/s", label="")
    # p2 = Plots.plot(sol, vars=[emf.i], ylabel="Current in A", label="")
    # Plots.plot(p1, p2, layout=(2,1))
    # Plots.savefig("dc_motor.png")
end

@testset "DC motor with speed sensor" begin
    R = 0.5
    L = 4.5e-3
    k = 0.5
    J = 0.02
    f = 0.01
    V_step = 10
    tau_L_step = -3
    @named ground = Ground()
    @named source = Voltage()
    @named voltage_step = Blocks.Step(height = V_step, start_time = 0, output__unit = u"V")
    @named R1 = Resistor(R = R)
    @named L1 = Inductor(L = L, i = 0.0)
    @named emf = EMF(k = k)
    @named fixed = Fixed()
    @named load = Torque(use_support = false)
    @named load_step = Blocks.Step(height = tau_L_step, start_time = 3, output__unit = u"N*m")
    @named inertia = Inertia(J = J)
    @named friction = Damper(d = f)
    @named speed_sensor = SpeedSensor()

    connections = [connect(fixed.flange, emf.support, friction.flange_b)
        connect(emf.flange, friction.flange_a, inertia.flange_a)
        connect(inertia.flange_b, load.flange)
        connect(inertia.flange_b, speed_sensor.flange)
        connect(load_step.output, load.tau)
        connect(voltage_step.output, source.V)
        connect(source.p, R1.p)
        connect(R1.n, L1.p)
        connect(L1.n, emf.p)
        connect(emf.n, source.n, ground.g)]

    @named model = ODESystem(connections, t,
        systems = [
            ground,
            voltage_step,
            source,
            R1,
            L1,
            emf,
            fixed,
            load,
            load_step,
            inertia,
            friction,
            speed_sensor,
        ])
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[], (0, 6.0))
    sol = solve(prob, Rodas4())

    @test sol.retcode == Success
    # EMF equations
    @test -0.5 .* sol[emf.i] == sol[emf.flange.tau]
    @test sol[emf.v] == 0.5 .* sol[emf.w]

    # test steady-state values
    dc_gain = [f/(k^2 + f * R) k/(k^2 + f * R); k/(k^2 + f * R) -R/(k^2 + f * R)]
    idx_t = findfirst(sol.t .> 2.5)
    @test sol[inertia.w][idx_t]≈(dc_gain * [V_step; 0])[2] rtol=1e-3
    @test sol[emf.i][idx_t]≈(dc_gain * [V_step; 0])[1] rtol=1e-3
    idx_t = findfirst(sol.t .> 5.5)
    @test sol[inertia.w][idx_t]≈(dc_gain * [V_step; -tau_L_step])[2] rtol=1e-3
    @test sol[emf.i][idx_t]≈(dc_gain * [V_step; -tau_L_step])[1] rtol=1e-3

    @test all(sol[inertia.w] .== sol[speed_sensor.w.u])

    prob = DAEProblem(sys, D.(unknowns(sys)) .=> 0.0, Pair[], (0, 6.0))
    sol = solve(prob, DFBDF())

    @test sol.retcode == Success
    # EMF equations
    @test -0.5 .* sol[emf.i] == sol[emf.flange.tau]
    @test sol[emf.v] == 0.5 .* sol[emf.w]
    # test steady-state values
    dc_gain = [f/(k^2 + f * R) k/(k^2 + f * R); k/(k^2 + f * R) -R/(k^2 + f * R)]
    idx_t = findfirst(sol.t .> 2.5)
    @test sol[inertia.w][idx_t]≈(dc_gain * [V_step; 0])[2] rtol=1e-3
    @test sol[emf.i][idx_t]≈(dc_gain * [V_step; 0])[1] rtol=1e-3
    idx_t = findfirst(sol.t .> 5.5)
    @test sol[inertia.w][idx_t]≈(dc_gain * [V_step; -tau_L_step])[2] rtol=1e-3
    @test sol[emf.i][idx_t]≈(dc_gain * [V_step; -tau_L_step])[1] rtol=1e-3
    #
    @test all(sol[inertia.w] .== sol[speed_sensor.w.u])
end

@testset "El. Heating Circuit" begin
    @named ground = Ground()
    @named source = Voltage()
    @named voltage_sine = Blocks.Sine(amplitude = 220, frequency = 1, output__unit = u"V")
    @named heating_resistor = HeatingResistor(R_ref = 100, alpha = 1e-3, T_ref = 293.15)
    @named thermal_conductor = ThermalConductor(G = 50)
    @named env = FixedTemperature(T = 273.15 + 20)
    connections = [connect(source.n, ground.g, heating_resistor.n)
        connect(source.p, heating_resistor.p)
        connect(voltage_sine.output, source.V)
        connect(heating_resistor.heat_port, thermal_conductor.port_a)
        connect(thermal_conductor.port_b, env.port)]

    @named model = ODESystem(connections, t,
        systems = [
            ground,
            voltage_sine,
            source,
            heating_resistor,
            thermal_conductor,
            env,
        ])
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[], (0, 6.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    ### TODO
    @test sol[source.v * source.i] == -sol[env.port.Q_flow]
end
