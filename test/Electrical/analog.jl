using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test

# using Plots

@parameters t

@testset "sensors" begin
    @named source = ConstantVoltage(V=10)
    @named resistor = Resistor(R=1)
    @named capacitor = Capacitor(C=1)
    @named ground = Ground()

    @named voltage_sensor = VoltageSensor()
    @named current_sensor = CurrentSensor()
    @named power_sensor = PowerSensor()

    connections = [
        connect(source.p, resistor.p)
        connect(resistor.n, current_sensor.p)
        connect(current_sensor.n, power_sensor.pc)
        connect(power_sensor.nc, capacitor.p)
        connect(capacitor.n, source.n, ground.g)
        connect(capacitor.p, voltage_sensor.p)
        connect(capacitor.n, voltage_sensor.n)
        connect(capacitor.p, power_sensor.pv)
        connect(capacitor.n, power_sensor.nv)
    ]

    @named model = ODESystem(connections, t; systems=[resistor, capacitor, source, ground, voltage_sensor, current_sensor, power_sensor])
    sys = structural_simplify(model)
    prob = ODAEProblem(sys, Pair[], (0.0, 10.0))
    @test_nowarn sol = solve(prob, Tsit5())

    # Plots.plot(sol; vars=[capacitor.v, voltage_sensor.v])
    # Plots.plot(sol; vars=[power_sensor.power, capacitor.i * capacitor.v])
    # Plots.plot(sol; vars=[resistor.i, current_sensor.i])
    @test sol[capacitor.v] ≈ sol[voltage_sensor.v] atol=1e-3
    @test sol[power_sensor.power] ≈ sol[capacitor.i * capacitor.v] atol=1e-3
    @test sol[resistor.i] ≈ sol[current_sensor.i] atol=1e-3
end

# simple voltage divider
@testset "voltage divider" begin
    @named source = ConstantVoltage(V=10)
    @named R1 = Resistor(R=1e3)
    @named R2 = Resistor(R=1e3)
    @named ground = Ground()

    connections = [
        connect(source.p, R1.p)
        connect(R1.n, R2.p)
        connect(R2.n, source.n, ground.g)
    ]

    @named model = ODESystem(connections, t, systems=[R1, R2, source, ground])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0, 2.0))
    sol = solve(prob, Rodas4()) # has no state; does not work with Tsit5
    @test sol[R1.p.v][end] ≈ 10 atol=1e-3
    @test sol[R1.n.v][end] ≈ 5 atol=1e-3
    @test sol[R2.n.v][end] ≈ 0 atol=1e-3
end

# simple RC
@testset "RC" begin
    @named source = ConstantVoltage(V=10)
    @named resistor = Resistor(R=1)
    @named capacitor = Capacitor(C=1)
    @named ground = Ground()

    connections = [
        connect(source.p, resistor.p)
        connect(resistor.n, capacitor.p)
        connect(capacitor.n, source.n, ground.g)
    ]

    @named model = ODESystem(connections, t; systems=[resistor, capacitor, source, ground])
    sys = structural_simplify(model)
    prob = ODAEProblem(sys, [capacitor.v => 0.0], (0.0, 10.0))
    sol = solve(prob, Tsit5())

    # Plots.plot(sol; vars=[source.v, capacitor.v])
    @test sol[capacitor.v][end] ≈ 10 atol=1e-3
end

# simple RL
@testset "RL" begin
    @named source = ConstantVoltage(V=10)
    @named resistor = Resistor(R=1)
    @named inductor = Inductor(L=1.0)
    @named ground = Ground()

    connections = [
        connect(source.p, resistor.p)
        connect(resistor.n, inductor.p)
        connect(inductor.n, source.n, ground.g)
    ]

    @named model = ODESystem(connections, t; systems=[resistor, inductor, source, ground])
    sys = structural_simplify(model)
    prob = ODAEProblem(sys, [inductor.i => 0.0], (0.0, 10.0))
    sol = solve(prob, Tsit5())

    # Plots.plot(sol; vars=[inductor.i, inductor.i])
    @test sol[inductor.i][end] ≈ 10 atol=1e-3
end

# RC with different voltage sources
@testset "RC with voltage sources" begin
    @named source_const = ConstantVoltage(V=10)
    @named source_sin = SineVoltage(offset=1, amplitude=10, frequency=2, start_time=0.5, phase=0)
    @named source_step = StepVoltage(offset=1, height=10, start_time=0.5)
    @named source_tri = TriangularVoltage(offset=1, start_time=0.5, amplitude=10, frequency=2)
    @named source_dsin = ExpSineVoltage(offset=1, amplitude=10, frequency=2, start_time=0.5, phase=0, damping=0.5)
    @named source_ramp = RampVoltage(offset=1, height=10, start_time=0.5, end_time=1.5)
    sources = [source_const, source_sin, source_step, source_tri, source_dsin, source_ramp]

    @named resistor = Resistor(R=1)
    @named capacitor = Capacitor(C=1)
    @named ground = Ground()

    for source in sources
        connections = [
            connect(source.p, resistor.p)
            connect(resistor.n, capacitor.p)
            connect(capacitor.n, source.n, ground.g)
        ]

        @named model = ODESystem(connections, t; systems=[resistor, capacitor, source, ground])
        sys = structural_simplify(model)
        prob = ODAEProblem(sys, [capacitor.v => 0.0], (0.0, 10.0))
        @test_nowarn sol = solve(prob, Tsit5())

        # Plots.plot(sol; vars=[source.v, capacitor.v])
    end
end

# RL with different voltage sources
@testset "RL with voltage sources" begin
    @named source_const = ConstantVoltage(V=10)
    @named source_sin = SineVoltage(offset=1, amplitude=10, frequency=2, start_time=0.5, phase=0)
    @named source_step = StepVoltage(offset=1, height=10, start_time=0.5)
    @named source_tri = TriangularVoltage(offset=1, start_time=0.5, amplitude=10, frequency=2)
    @named source_dsin = ExpSineVoltage(offset=1, amplitude=10, frequency=2, start_time=0.5, phase=0, damping=0.5)
    @named source_ramp = RampVoltage(offset=1, height=10, start_time=0.5, end_time=1.5)
    sources = [source_const, source_sin, source_step, source_tri, source_dsin, source_ramp]

    @named resistor = Resistor(R=1.0)
    @named inductor = Inductor(L=1.0)
    @named ground = Ground()

    for source in sources
        connections = [
            connect(source.p, resistor.p)
            connect(resistor.n, inductor.p)
            connect(inductor.n, source.n, ground.g)
        ]

        @named model = ODESystem(connections, t; systems=[resistor, inductor, source, ground])
        sys = structural_simplify(model)
        prob = ODAEProblem(sys, [inductor.i => 0.0], (0.0, 10.0))
        @test_nowarn sol = solve(prob, Tsit5())

        # Plots.plot(sol; vars=[source.i, inductor.i])
    end
end

# RC with different current sources
@testset "RC with current sources" begin
    @named source_const = ConstantCurrent(I=10)
    @named source_sin = SineCurrent(offset=1, amplitude=10, frequency=2, start_time=0.5, phase=0)
    @named source_step = StepCurrent(offset=1, height=10, start_time=0.5)
    @named source_tri = TriangularCurrent(offset=1, start_time=0.5, amplitude=10, frequency=2)
    @named source_dsin = ExpSineCurrent(offset=1, amplitude=10, frequency=2, start_time=0.5, phase=0, damping=0.5)
    @named source_ramp = RampCurrent(offset=1, height=10, start_time=0.5, end_time=1.5)
    sources = [source_const, source_sin, source_step, source_tri, source_dsin, source_ramp]

    @named resistor = Resistor(R=1)
    @named capacitor = Capacitor(C=1)
    @named ground = Ground()

    for source in sources
        connections = [
            connect(source.p, resistor.p)
            connect(resistor.n, capacitor.p)
            connect(capacitor.n, source.n, ground.g)
        ]

        @named model = ODESystem(connections, t; systems=[resistor, capacitor, source, ground])
        sys = structural_simplify(model)
        prob = ODAEProblem(sys, [capacitor.v => 0.0], (0.0, 10.0))
        @test_nowarn sol = solve(prob, Tsit5())

        # Plots.plot(sol; vars=[source.v, capacitor.v])
    end
end
@testset "Integrator" begin
    R=1e3
    f=1
    Vin=5
    @named ground = Ground()
    @named R1 = Resistor(R=R)
    @named R2 = Resistor(R=100*R)
    @named C1 = Capacitor(C=1/(2 * pi * f * R))
    @named opamp = IdealOpAmp()
    @named square = SquareVoltage(amplitude=Vin)
    @named sensor = VoltageSensor()
    
    connections = [
        connect(square.p, R1.p)
        connect(R1.n, C1.p, R2.p, opamp.n1)
        connect(opamp.p2, C1.n, R2.n)
        connect(opamp.p1, ground.g, opamp.n2, square.n)
        connect(opamp.p2, sensor.p)
        connect(sensor.n, ground.g)
    ]
    @named model = ODESystem(connections, t, systems = [R1, R2, opamp, square, C1, ground, sensor])
    sys = structural_simplify(model)
    u0 = [
        C1.v => 0.0
        R1.v => 0.0
    ]
    prob = ODEProblem(sys, u0, (0, 100.0))
    sol = solve(prob, Rodas4())
    @test sol[opamp.v2] == sol[-C1.v] # Not a great one however. Rely on the plot
    @test sol[opamp.p2.v] == sol[sensor.v] 

    # plot(sol, vars=[sensor.v, square.v, C1.v])
end
