using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkitStandardLibrary.Blocks: Step, Constant, Sine, Cosine, ExpSine
# using Plots

@parameters t

@testset "sensors" begin
    @named source    = Sine(offset=1, amplitude=10, frequency=5)
    @named voltage   = Voltage()
    @named resistor  = Resistor(R=1)
    @named capacitor = Capacitor(C=1)
    @named ground    = Ground()

    @named voltage_sensor = VoltageSensor()
    @named current_sensor = CurrentSensor()
    @named power_sensor   = PowerSensor()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, resistor.p)
        connect(resistor.n, current_sensor.p)
        connect(current_sensor.n, power_sensor.pc)
        connect(power_sensor.nc, capacitor.p)
        connect(capacitor.n, voltage.n, ground.g)
        connect(capacitor.p, voltage_sensor.p)
        connect(capacitor.n, voltage_sensor.n)
        connect(capacitor.p, power_sensor.pv)
        connect(capacitor.n, power_sensor.nv)
    ]

    @named model = ODESystem(connections, t; systems=[resistor, capacitor, source, voltage, ground, voltage_sensor, current_sensor, power_sensor])
    sys = structural_simplify(model)
    prob = ODAEProblem(sys, Pair[], (0.0, 10.0))
    sol = solve(prob, Tsit5())

    # Plots.plot(sol; vars=[capacitor.v, voltage_sensor.v])
    # Plots.plot(sol; vars=[power_sensor.power, capacitor.i * capacitor.v])
    # Plots.plot(sol; vars=[resistor.i, current_sensor.i])
    @test sol[capacitor.v] ≈ sol[voltage_sensor.v] atol=1e-3
    @test sol[power_sensor.power] ≈ sol[capacitor.i * capacitor.v] atol=1e-3
    @test sol[resistor.i] ≈ sol[current_sensor.i] atol=1e-3
end

# simple voltage divider
@testset "voltage divider" begin
    @named source  = Constant(k=10)
    @named voltage = Voltage()
    @named R1      = Resistor(R=1e3)
    @named R2      = Resistor(R=1e3)
    @named ground  = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, R1.p)
        connect(R1.n, R2.p)
        connect(R2.n, voltage.n, ground.g)
    ]

    @named model = ODESystem(connections, t, systems=[R1, R2, source,voltage, ground])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0, 2.0))
    sol = solve(prob, Rodas4()) # has no state; does not work with Tsit5
    @test sol[R1.p.v][end] ≈ 10 atol=1e-3
    @test sol[R1.n.v][end] ≈ 5 atol=1e-3
    @test sol[R2.n.v][end] ≈ 0 atol=1e-3
end

# simple RC
@testset "RC" begin
    @named source    = Constant(k=10)
    @named voltage   = Voltage()
    @named resistor  = Resistor(R=1)
    @named capacitor = Capacitor(C=1)
    @named ground    = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, resistor.p)
        connect(resistor.n, capacitor.p)
        connect(capacitor.n, voltage.n, ground.g)
    ]

    @named model = ODESystem(connections, t; systems=[resistor, capacitor, source, voltage, ground])
    sys = structural_simplify(model)
    prob = ODAEProblem(sys, [capacitor.v => 0.0], (0.0, 10.0))
    sol = solve(prob, Tsit5())

    # Plots.plot(sol; vars=[source.v, capacitor.v])
    @test sol[capacitor.v][end] ≈ 10 atol=1e-3
end

# simple RL
@testset "RL" begin
    @named source   = Constant(k=10)
    @named voltage  = Voltage()
    @named resistor = Resistor(R=1)
    @named inductor = Inductor(L=1.0)
    @named ground   = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, resistor.p)
        connect(resistor.n, inductor.p)
        connect(inductor.n, voltage.n, ground.g)
    ]

    @named model = ODESystem(connections, t; systems=[resistor, inductor, source, voltage, ground])
    sys = structural_simplify(model)
    prob = ODAEProblem(sys, [inductor.i => 0.0], (0.0, 10.0))
    sol = solve(prob, Tsit5())

    # Plots.plot(sol; vars=[inductor.i, inductor.i])
    @test sol[inductor.i][end] ≈ 10 atol=1e-3
end

# RC with current sources
@testset "RC with current sources" begin
    start_time = 2
    @named current   = Current()
    @named source    = Step(start_time=2)
    @named resistor  = Resistor(R=1)
    @named capacitor = Capacitor(C=1)
    @named ground    = Ground()

    connections = [
        connect(source.output, current.I)
        connect(current.p, resistor.n)
        connect(capacitor.n, resistor.p)
        connect(capacitor.p, current.n, ground.g)
    ]

    @named model = ODESystem(connections, t; systems=[ground, resistor,  current, capacitor, source])
    sys = structural_simplify(model)
    prob = ODAEProblem(sys, [capacitor.v => 0.0], (0.0, 10.0))
    sol = solve(prob, Tsit5())
    y(x, st) = (x .> st) .* abs.(collect(x) .- st)
    @test sum(reduce(vcat, sol.u) .- y(sol.t, start_time)) ≈ 0 atol=1e-4 

end
#=
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
        connect(R1.n, C1.n, R2.p, opamp.n1)
        connect(opamp.p2, C1.p, R2.n)
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
    @test sol[opamp.v2] == sol[C1.v] # Not a great one however. Rely on the plot
    @test sol[opamp.p2.v] == sol[sensor.v] 

    # plot(sol, vars=[sensor.v, square.v, C1.v])
end
=#
_step(x, h, st) = ifelse(x < st, 0, h)
_cos_wave(x, f, A, st, ϕ) = A*cos(2*π*f*(x-st) + ϕ)
# _ramp(x, st, d, h) = ifelse(x < st, 0,
#                                   ifelse(x < (st + d), (x - st) * h / d, h))
_sine_wave(x, f, A, st, ϕ) = A*sin(2*π*f*(x - st) + ϕ)
_damped_sine_wave(x, f, A, st, ϕ, d) = exp((st-x)*d)*A*sin(2*π*f*(x-st) + ϕ)

@testset "Voltage function generators" begin
    st, o, h, f, A, et, ϕ, d, δ = 0.7, 1.25, 3, 2, 2.5, 2, π/4, 0.1, 0.0001

    @named res     = Resistor(R=1)
    @named cap     = Capacitor(C=1)
    @named ground  = Ground()
    @named voltage = Voltage()
    @named voltage_sensor = VoltageSensor()
    @named step   = Step(start_time=st, offset=o, height=h)
    @named cosine = Cosine(offset=o, amplitude=A, frequency=f, start_time=st, phase=ϕ)
    @named sine   = Sine(offset=o, amplitude=A, frequency=f, start_time=st, phase=ϕ)
    @named damped_sine = ExpSine(offset=o, amplitude=A, frequency=f, start_time=st, phase=ϕ, damping=d)
    # @named vsquare = Square(offset=o, start_time=st, amplitude=A, frequency=f)
    # @named vtri = Triangular(offset=o, start_time=st, amplitude=A, frequency=f)
    # @named vsawtooth = SawTooth(amplitude=A, start_time=st, frequency=f, offset=o)
    # @named vramp   = Ramp(offset=o, start_time=st, duration=et-st, height=h)

    sources = [step, cosine, sine, damped_sine] #, vtri, vsquare, vdamped_sine]
    waveforms(i, x) = getindex([o .+  _step.(x, h, st),
                                o .+ (x .> st) .* _cos_wave.(x, f, A, st, ϕ),
                                o .+ (x .> st) .* _sine_wave.(x, f, A, st, ϕ),
                                o .+ (x .> st) .* _damped_sine_wave.(x, f, A, st, ϕ, d)], i)
                                # o .+ (x .> st) .* _triangular_wave.(x, δ, f, A, st),
                                # o .+ (x .> st) .* _square_wave.(x, δ, f, A, st),
                                # o .+ (x .> st). * _sawtooth_wave.(x, δ, f, A, st),
                                # o .+ _ramp.(x, st, et, h)], i)

    for i in 1:length(sources)
        source = sources[i]
        eqs = [
            connect(source.output, voltage.V)
            connect(voltage.p, voltage_sensor.p, res.p)
            connect(res.n, cap.p)
            connect(ground.g, voltage_sensor.n, voltage.n, cap.n)
        ]
        @named vmodel = ODESystem(eqs, t, systems = [voltage_sensor, res, cap, source, voltage, ground])
        vsys = structural_simplify(vmodel)

        u0 = []

        prob = ODAEProblem(vsys, u0, (0, 10.0))
        sol = solve(prob, dt=0.1, Tsit5())

        @test sol[voltage.V.u] ≈ waveforms(i, sol.t) atol=1e-1
        @test sol[voltage.p.v] ≈ sol[voltage.V.u]
        # For visual inspection
        # plt = plot(sol; vars=[vsource.v])
        # savefig(plt, "test_voltage_$(Symbolics.getname(vsource))")
    end
end

@testset "Current function generators" begin
    st, o, h, f, A, et, ϕ, d, δ = 0.7, 1.25, 3, 2, 2.5, 2, π/4, 0.1, 0.0001

    @named ground = Ground()
    @named res = Resistor(R=1.0)
    @named cap = Capacitor(C=1)
    @named current_sensor = CurrentSensor()
    @named current = Current()
    @named step = Step(start_time=st, offset=o, height=h)
    @named cosine = Cosine(offset=o, amplitude=A, frequency=f, start_time=st, phase=ϕ)
    @named sine = Sine(offset=o, amplitude=A, frequency=f, start_time=st, phase=ϕ)
    @named damped_sine = ExpSine(offset=o, amplitude=A, frequency=f, start_time=st, phase=ϕ, damping=d)
    # @named isquare = Square(offset=o, start_time=st, amplitude=A, frequency=f)
    # @named itri = Triangular(offset=o, start_time=st, amplitude=A, frequency=f)
    # @named isawtooth = SawTooth(amplitude=A, start_time=st, frequency=f, offset=o)
    # @named iramp = Ramp(offset=o, start_time=st, duration=et-st, height=h)

    sources = [step, cosine, sine, damped_sine] #, itri, isquare, idamped_sine]
    waveforms(i, x) = getindex([o .+  _step.(x, h, st),
                                o .+ (x .> st) .* _cos_wave.(x, f, A, st, ϕ),
                                o .+ (x .> st) .* _sine_wave.(x, f, A, st, ϕ),
                                o .+ (x .> st) .* _damped_sine_wave.(x, f, A, st, ϕ, d)], i)
                                # o .+ (x .> st) .* _triangular_wave.(x, δ, f, A, st),
                                # o .+ (x .> st) .* _square_wave.(x, δ, f, A, st),
                                # # o .+ (x .> st). * _sawtooth_wave.(x, δ, f, A, st),
                                # o .+ _ramp.(x, δ, st, et, h)], i)

    for i in 1:length(sources)
        source = sources[i]
        eqs = [
            connect(source.output, current.I)
            connect(current.p, current_sensor.n)
            connect(current_sensor.p, res.p)
            connect(res.n, cap.p)
            connect(current.n, ground.g, cap.n)
        ]
        @named model = ODESystem(eqs, t, systems = [current_sensor, source, current, res, cap, ground])
        isys = structural_simplify(model)

        u0 = []

        prob = ODAEProblem(isys, u0, (0, 10.0))
        sol = solve(prob, Tsit5())

        @test sol[current.I.u] ≈ waveforms(i, sol.t) atol=1e-1
        @test sol[current.I.u] ≈ sol[current.p.i] atol=1e-1
        # For visual inspection
        # plt = plot(sol)
        # savefig(plt, "test_current_$(Symbolics.getname(isource))")
    end
end