using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Blocks: Step,
                                             Constant, Sine, Cosine, ExpSine, Ramp,
                                             Square, Triangular
using ModelingToolkitStandardLibrary.Blocks: square, triangular
using ModelingToolkitStandardLibrary.Thermal: FixedTemperature
using OrdinaryDiffEq: ReturnCode.Success

# using Plots

@testset "sensors" begin
    @named source = Sine(offset = 1, amplitude = 10, frequency = 5)
    @named voltage = Voltage()
    @named resistor = Resistor(R = 1)
    @named capacitor = Capacitor(C = 1, v = 0.0)
    @named ground = Ground()

    @named voltage_sensor = VoltageSensor()
    @named current_sensor = CurrentSensor()
    @named power_sensor = PowerSensor()

    connections = [connect(source.output, voltage.V)
                   connect(voltage.p, resistor.p)
                   connect(resistor.n, current_sensor.p)
                   connect(current_sensor.n, power_sensor.pc)
                   connect(power_sensor.nc, capacitor.p)
                   connect(capacitor.n, voltage.n, ground.g)
                   connect(capacitor.p, voltage_sensor.p)
                   connect(capacitor.n, voltage_sensor.n)
                   connect(capacitor.p, power_sensor.pv)
                   connect(capacitor.n, power_sensor.nv)]

    @named model = ODESystem(connections, t;
        systems = [
            resistor,
            capacitor,
            source,
            voltage,
            ground,
            voltage_sensor,
            current_sensor,
            power_sensor
        ])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, [], (0.0, 10.0))
    sol = solve(prob, Tsit5())

    # Plots.plot(sol; vars=[capacitor.v, voltage_sensor.v])
    # Plots.plot(sol; vars=[power_sensor.power, capacitor.i * capacitor.v])
    # Plots.plot(sol; vars=[resistor.i, current_sensor.i])
    @test SciMLBase.successful_retcode(sol)
    @test sol[capacitor.v]≈sol[voltage_sensor.v] atol=1e-3
    @test sol[power_sensor.power]≈sol[capacitor.i * capacitor.v] atol=1e-3
    @test sol[resistor.i]≈sol[current_sensor.i] atol=1e-3
end

# simple voltage divider
@testset "voltage divider with a short branch" begin
    @named source = Constant(k = 10)
    @named voltage = Voltage()
    @named R0 = Resistor(R = 1e3)
    @named R1 = Resistor(R = 1e3)
    @named R2 = Resistor(R = 1e3)
    @named ground = Ground()
    @named short = Short()

    connections = [connect(source.output, voltage.V)
                   connect(voltage.p, R1.p)
                   connect(R1.n, short.p, R0.p)
                   connect(short.n, R0.n, R2.p)
                   connect(R2.n, voltage.n, ground.g)]

    @named model = ODESystem(connections, t,
        systems = [R0, R1, R2, source, short, voltage, ground]; guesses = [R2.v => 0.0])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, [], (0, 2.0))
    sol = solve(prob, Rodas4()) # has no state; does not work with Tsit5
    @test SciMLBase.successful_retcode(sol)
    @test sol[short.v] == sol[R0.v] == zeros(length(sol.t))
    @test sol[R0.i] == zeros(length(sol.t))
    @test sol[R1.p.v][end]≈10 atol=1e-3
    @test sol[R1.n.v][end]≈5 atol=1e-3
    @test sol[R2.n.v][end]≈0 atol=1e-3
end

# simple RC
@testset "RC" begin
    @named source = Constant(k = 10)
    @named voltage = Voltage()
    @named resistor = Resistor(R = 1)
    @named capacitor = Capacitor(C = 1, v = 0.0)
    @named ground = Ground()

    connections = [connect(source.output, voltage.V)
                   connect(voltage.p, resistor.p)
                   connect(resistor.n, capacitor.p)
                   connect(capacitor.n, voltage.n, ground.g)]

    @named model = ODESystem(connections, t;
        systems = [resistor, capacitor, source, voltage, ground])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0.0, 10.0))
    sol = solve(prob, Tsit5())

    # Plots.plot(sol; vars=[source.v, capacitor.v])
    @test SciMLBase.successful_retcode(sol)
    @test sol[capacitor.v][end]≈10 atol=1e-3
end

# simple RL
@testset "RL" begin
    @named source = Constant(k = 10)
    @named voltage = Voltage()
    @named resistor = Resistor(R = 1)
    @named inductor = Inductor(L = 1.0, i = 0.0)
    @named ground = Ground()

    connections = [connect(source.output, voltage.V)
                   connect(voltage.p, resistor.p)
                   connect(resistor.n, inductor.p)
                   connect(inductor.n, voltage.n, ground.g)]

    @named model = ODESystem(connections, t;
        systems = [resistor, inductor, source, voltage, ground])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0.0, 10.0))
    sol = solve(prob, Tsit5())

    # Plots.plot(sol; vars=[inductor.i, inductor.i])
    @test SciMLBase.successful_retcode(sol)
    @test sol[inductor.i][end]≈10 atol=1e-3
end

@testset "RC with voltage sources" begin
    R, C = 1, 1
    @named voltage = Voltage()
    @named source_const = Constant(k = 10)
    @named source_sin = Sine(offset = 1, amplitude = 10, frequency = 2, start_time = 0.5,
        phase = 0)
    @named source_step = Step(offset = 1, height = 10, start_time = 0.5)
    @named source_tri = Triangular(offset = 1, start_time = 0.5, amplitude = 10,
        frequency = 2)
    @named source_dsin = ExpSine(offset = 1, amplitude = 10, frequency = 2,
        start_time = 0.5, phase = 0, damping = 0.5)
    @named source_ramp = Ramp(offset = 1, height = 10, start_time = 0.5, duration = 1)
    sources = [source_const, source_sin, source_step, source_tri, source_dsin, source_ramp]

    @named resistor = Resistor(; R)
    @named capacitor = Capacitor(; C, v = 0.0)
    @named ground = Ground()

    for source in sources
        connections = [connect(source.output, voltage.V)
                       connect(voltage.p, resistor.p)
                       connect(resistor.n, capacitor.p)
                       connect(capacitor.n, voltage.n, ground.g)]

        @named model = ODESystem(connections, t;
            systems = [resistor, capacitor, source, ground, voltage])
        sys = structural_simplify(model)
        prob = ODEProblem(sys, Pair[], (0.0, 10.0))
        sol = solve(prob, Tsit5())
        @test SciMLBase.successful_retcode(sol)
        sol = solve(prob, Rodas4())
        @test SciMLBase.successful_retcode(sol)

        # Plots.plot(sol; vars=[voltage.v, capacitor.v])
    end
end

# RC with current sources
@testset "RC with current sources" begin
    start_time = 2
    @named current = Current()
    @named source = Step(start_time = 2)
    @named resistor = Resistor(R = 1)
    @named capacitor = Capacitor(C = 1, v = 0.0)
    @named ground = Ground()

    connections = [connect(source.output, current.I)
                   connect(current.p, resistor.n)
                   connect(capacitor.n, resistor.p)
                   connect(capacitor.p, current.n, ground.g)]

    @named model = ODESystem(connections, t;
        systems = [ground, resistor, current, capacitor, source])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0.0, 10.0))
    sol = solve(prob, Tsit5())
    y(x, st) = (x .> st) .* abs.(collect(x) .- st)
    @test SciMLBase.successful_retcode(sol)
    @test sum(reduce(vcat, sol[capacitor.v]) .- y(sol.t, start_time))≈0 atol=1e-2
end

@testset "Integrator" begin
    R = 1e3
    f = 1
    Vin = 5
    @named ground = Ground()
    @named R1 = Resistor(R = R)
    @named R2 = Resistor(R = 100 * R)
    @named C1 = Capacitor(C = 1 / (2 * pi * f * R), v = 0.0)
    @named opamp = IdealOpAmp()
    @named square_source = Square(amplitude = Vin)
    @named voltage = Voltage()
    @named sensor = VoltageSensor()

    connections = [connect(square_source.output, voltage.V)
                   connect(voltage.p, R1.p)
                   connect(R1.n, C1.n, R2.p, opamp.n1)
                   connect(opamp.p2, C1.p, R2.n)
                   connect(opamp.p1, ground.g, opamp.n2, voltage.n)
                   connect(opamp.p2, sensor.p)
                   connect(sensor.n, ground.g)]
    @named model = ODESystem(connections, t,
        systems = [
            R1,
            R2,
            opamp,
            square_source,
            voltage,
            C1,
            ground,
            sensor
        ])
    sys = structural_simplify(model)
    u0 = [C1.v => 0.0
          R1.v => 0.0]
    prob = ODEProblem(sys, u0, (0, 100.0))
    sol = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)
    @test sol[opamp.v2] == sol[C1.v] # Not a great one however. Rely on the plot
    @test sol[opamp.p2.v] == sol[sensor.v]

    # plot(sol, vars=[sensor.v, square.v, C1.v])
end

_step(x, h, st) = ifelse(x < st, 0, h)
_cos_wave(x, f, A, st, ϕ) = A * cos(2 * π * f * (x - st) + ϕ)
_ramp(x, st, d, h) = ifelse(x < st, 0,
    ifelse(x < (st + d), (x - st) * h / d, h))
_sine_wave(x, f, A, st, ϕ) = A * sin(2 * π * f * (x - st) + ϕ)
_damped_sine_wave(x, f, A, st, ϕ, d) = exp((st - x) * d) * A * sin(2 * π * f * (x - st) + ϕ)

@testset "Voltage function generators" begin
    st, o, h, f, A, et, ϕ, d, δ = 0.7, 1.25, 3, 2, 2.5, 2, π / 4, 0.1, 0.0001

    @named res = Resistor(R = 1)
    @named cap = Capacitor(C = 1, v = 0.0)
    @named ground = Ground()
    @named voltage = Voltage()
    @named voltage_sensor = VoltageSensor()
    @named step = Step(start_time = st, offset = o, height = h)
    @named cosine = Cosine(offset = o, amplitude = A, frequency = f, start_time = st,
        phase = ϕ)
    @named sine = Sine(offset = o, amplitude = A, frequency = f, start_time = st, phase = ϕ)
    @named damped_sine = ExpSine(offset = o, amplitude = A, frequency = f, start_time = st,
        phase = ϕ, damping = d)
    @named ramp = Ramp(offset = o, start_time = st, duration = et - st, height = h)
    @named vsquare = Square(offset = o, start_time = st, amplitude = A, frequency = f)
    @named tri = Triangular(offset = o, start_time = st, amplitude = A, frequency = f)
    # @named vsawtooth = SawTooth(amplitude=A, start_time=st, frequency=f, offset=o)

    sources = [step, cosine, sine, damped_sine, ramp, tri, vsquare] #, vsawtooth]
    function waveforms(i, x)
        getindex(
            [o .+ _step.(x, h, st),
                o .+ (x .> st) .* _cos_wave.(x, f, A, st, ϕ),
                o .+ (x .> st) .* _sine_wave.(x, f, A, st, ϕ),
                o .+ (x .> st) .* _damped_sine_wave.(x, f, A, st, ϕ, d),
                o .+ _ramp.(x, st, (et - st), h),
                triangular.(x, f, A, o, st),
                square.(x, f, A, o, st)],
            i)
    end
    # o .+ (x .> st). * _sawtooth_wave.(x, δ, f, A, st),

    for i in 1:lastindex(sources)
        source = sources[i]
        @info "Testing Voltage with $(nameof(source)) source"
        eqs = [connect(source.output, voltage.V)
               connect(voltage.p, voltage_sensor.p, res.p)
               connect(res.n, cap.p)
               connect(ground.g, voltage_sensor.n, voltage.n, cap.n)]
        @named vmodel = ODESystem(eqs, t,
            systems = [
                voltage_sensor,
                res,
                cap,
                source,
                voltage,
                ground
            ])
        vsys = structural_simplify(vmodel)

        u0 = [cap.v => 0.0]

        prob = ODEProblem(vsys, u0, (0, 10.0))
        sol = solve(prob, dt = 0.1, Tsit5())

        @test SciMLBase.successful_retcode(sol)
        @test sol[voltage.V.u]≈waveforms(i, sol.t) atol=1e-1
        @test sol[voltage.p.v] ≈ sol[voltage.V.u]
        # For visual inspection
        # plt = plot(sol; vars=[voltage.v])
        # savefig(plt, "test_voltage_$(source.name)")
    end
end

@testset "Current function generators" begin
    st, o, h, f, A, et, ϕ, d, δ = 0.7, 1.25, 3, 2, 2.5, 2, π / 4, 0.1, 0.0001

    @named ground = Ground()
    @named res = Resistor(R = 1.0)
    @named cap = Capacitor(C = 1, v = 0.0)
    @named current_sensor = CurrentSensor()
    @named current = Current()
    @named step = Step(start_time = st, offset = o, height = h)
    @named cosine = Cosine(offset = o, amplitude = A, frequency = f, start_time = st,
        phase = ϕ)
    @named sine = Sine(offset = o, amplitude = A, frequency = f, start_time = st, phase = ϕ)
    @named damped_sine = ExpSine(offset = o, amplitude = A, frequency = f, start_time = st,
        phase = ϕ, damping = d)
    @named ramp = Ramp(offset = o, start_time = st, duration = et - st, height = h)
    @named vsquare = Square(offset = o, start_time = st, amplitude = A, frequency = f)
    @named tri = Triangular(offset = o, start_time = st, amplitude = A, frequency = f)
    # @named isawtooth = SawTooth(amplitude=A, start_time=st, frequency=f, offset=o)

    sources = [step, cosine, sine, damped_sine, ramp, tri, vsquare] #, idamped_sine]
    function waveforms(i, x)
        getindex(
            [o .+ _step.(x, h, st),
                o .+ (x .> st) .* _cos_wave.(x, f, A, st, ϕ),
                o .+ (x .> st) .* _sine_wave.(x, f, A, st, ϕ),
                o .+ (x .> st) .* _damped_sine_wave.(x, f, A, st, ϕ, d),
                o .+ _ramp.(x, st, (et - st), h),
                triangular.(x, f, A, o, st),
                square.(x, f, A, o, st)],
            i)
    end
    # # o .+ (x .> st). * _sawtooth_wave.(x, δ, f, A, st)

    for i in 1:lastindex(sources)
        source = sources[i]
        @info "Testing Current with $(nameof(source)) source"
        eqs = [connect(source.output, current.I)
               connect(current.p, current_sensor.n)
               connect(current_sensor.p, res.p)
               connect(res.n, cap.p)
               connect(current.n, ground.g, cap.n)]
        @named model = ODESystem(eqs, t,
            systems = [
                current_sensor,
                source,
                current,
                res,
                cap,
                ground
            ])
        isys = structural_simplify(model)

        u0 = [cap.v => 0.0]

        prob = ODEProblem(isys, u0, (0, 10.0))
        sol = solve(prob, dt = 0.1, Tsit5())

        @test SciMLBase.successful_retcode(sol)
        @test sol[current.I.u]≈waveforms(i, sol.t) atol=1e-1
        @test sol[current.I.u]≈sol[current.p.i] atol=1e-1
        # For visual inspection
        # plt = plot(sol)
        # savefig(plt, "test_current_$(source.name)")
    end
end

@testset "Diode component test" begin
    # Parameter values 
    R = 1.0
    C = 1.0
    V = 10.0
    n = 1.0
    Is = 1e-3
    f = 1.0

    # Components
    @named resistor = Resistor(R = R)
    @named capacitor = Capacitor(C = C, v = 0.0)
    @named source = Voltage()
    @named diode = Diode(n = n, Is = Is)
    @named ac = Sine(frequency = f, amplitude = V)
    @named ground = Ground()

    # Connections
    connections = [connect(ac.output, source.V)
                   connect(source.p, diode.p)
                   connect(diode.n, resistor.p)
                   connect(resistor.n, capacitor.p)
                   connect(capacitor.n, source.n, ground.g)]

    # Model
    @named model = ODESystem(connections, t;
        systems = [resistor, capacitor, source, diode, ac, ground])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0.0, 10.0))
    sol = solve(prob)

    # Extract solutions for testing
    diode_voltage = sol[diode.v]
    diode_current = sol[diode.i]
    resistor_current = sol[resistor.i]
    capacitor_voltage = sol[capacitor.v]

    # Tests
    @test all(diode_current .>= -Is)
    @test capacitor_voltage[end].≈V rtol=3e-1

    # For visual inspection
    # plt = plot(sol; idxs = [diode.i, resistor.i, capacitor.v],
    #     size = (800, 600), dpi = 300,
    #     labels = ["Diode Current" "Resistor Current" "Capacitor Voltage"],
    #     title = "Diode Test")
    # savefig(plt, "diode_test")
end

@testset "Diode with temperature dependency component test" begin
    # Parameter values
    R = 1.0
    C = 1.0
    V = 10.0
    T = 300.0 # Ambient temperature in Kelvin
    n = 2.0
    Is = 1e-6
    f = 1.0

    # Components
    @named resistor = Resistor(R = R)
    @named capacitor = Capacitor(C = C, v = 0.0)
    @named source = Voltage()
    @named heating_diode = Diode(n = n, Is = Is, T_dep = true)
    @named ac = Sine(frequency = f, amplitude = V)
    @named ground = Ground()
    @named temp = FixedTemperature(T = T)

    # Connections
    connections = [connect(ac.output, source.V),
        connect(source.p, heating_diode.p),
        connect(heating_diode.n, resistor.p),
        connect(resistor.n, capacitor.p),
        connect(capacitor.n, ground.g),
        connect(source.n, ground.g),
        connect(temp.port, heating_diode.port)]

    # Model
    @named model = ODESystem(connections, t;
        systems = [resistor, capacitor, source, heating_diode, ac, ground, temp])
    sys = structural_simplify(model)
    prob = ODEProblem(sys, Pair[], (0.0, 10.0))
    sol = solve(prob)

    # Extract solutions for testing
    diode_voltage = sol[heating_diode.v]
    diode_current = sol[heating_diode.i]
    resistor_current = sol[resistor.i]
    capacitor_voltage = sol[capacitor.v]

    # Expected thermal voltage at given temperature
    k = 1.380649e-23  # Boltzmann constant (J/K)
    q = 1.602176634e-19  # Elementary charge (C)

    # Tests
    @test all(diode_current .>= -Is) # Diode current should not exceed reverse saturation
    @test capacitor_voltage[end]≈V rtol=3e-1 # Final capacitor voltage close to input voltage

    # For visual inspection
    # plt = plot(sol; idxs = [heating_diode.i, resistor.i, capacitor.v],
    #     size = (800, 600), dpi = 300,
    #     labels = ["Diode Current" "Resistor Current" "Capacitor Voltage"],
    #     title = "Diode Test")
    # savefig(plt, "heating_diode_test")

    # Remake model with higher amb. temperature, final capacitor voltage should be lower
    T = 400.0
    model = remake(prob; p = [temp.T => T])
    sol = solve(model)
    @test SciMLBase.successful_retcode(sol)
    @test sol[capacitor.v][end] < capacitor_voltage[end]
end

@testset "VariableResistor with Temperature Dependency" begin
    R_ref = 2.0
    R_const = 1.0

    # Define the RC model as described
    @mtkmodel RC begin
        @parameters begin
            R = R_ref  # Variable resistance reference value
            C = 1.0   # Capacitance
            k = 10.0  # Voltage source scaling factor
            f = 0.2   # Frequency of sine input
            T = 300.0 # Ambient temperature in Kelvin
        end
        @components begin
            res_input = Sine(frequency = f, amplitude = 1.0, offset = 0.0)
            volt_input = Constant(k = 1.0)
            resistor = VariableResistor(R_ref = R_ref, R_const = R_const, T_dep = true)
            capacitor = Capacitor(C = C, v = 0.0)
            source = Voltage()
            temp = FixedTemperature(T = T)
            ground = Ground()
        end
        @equations begin
            connect(temp.port, resistor.port)
            connect(res_input.output, resistor.position)
            connect(volt_input.output, source.V)
            connect(source.p, resistor.p)
            connect(resistor.n, capacitor.p)
            connect(capacitor.n, source.n, ground.g)
        end
    end

    # Build and solve the system
    @mtkbuild sys = RC()
    prob = ODEProblem(sys, [], (0.0, 10.0); guesses = [sys.resistor.i => 0.0]) # No state variables initially
    sol = solve(prob)

    # Perform Tests
    resistor_resistance = sol[sys.resistor.R]
    capacitor_voltage = sol[sys.capacitor.v]

    @test SciMLBase.successful_retcode(sol) # Ensure the simulation is successful
    @test all(resistor_resistance .>= R_const) # Resistance should be >= constant value
    @test maximum(resistor_resistance) ≤ R_const + R_ref # Maximum resistance when pos=1 (R_const + R_ref)
    @test all(capacitor_voltage .>= 0.0) # Capacitor voltage should not be negative

    # For visual inspection
    # plt = plot(sol; vars = [sys.resistor.R, sys.capacitor.v],
    #     size = (800, 600), dpi = 300,
    #     labels = ["Variable Resistor Resistance" "Capacitor Voltage"],
    #     title = "RC Circuit Test with VariableResistor")
    # savefig(plt, "rc_circuit_test_variable_resistor")
end
