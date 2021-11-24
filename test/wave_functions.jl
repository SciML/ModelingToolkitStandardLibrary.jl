using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Blocks: _step, _square_wave, _triangular_wave,
                                 _cos_wave, _damped_sine_wave, _ramp

@parameters t
@named ground   = Ground()
@named resistor = Resistor()

st, o, h, f, A, et, ϕ, d, δ = 0.7, 1.25, 3, 2, 2.5, 2.5, π/4, 0.1, 0.0001
waveforms(i, x) = getindex([o .+ (x .> st) .* _triangular_wave.(x, δ, f, A, st),
                            o .+ (x .> st) .* _square_wave.(x, δ, f, A, st),
                            o .+  _step.(x, δ, h, st),
                            # o .+ (x .> st). * _sawtooth_wave.(x, δ, f, A, st),
                            o .+ (x .> st) .* _cos_wave.(x, f, A, st, ϕ),
                            o .+ (x .> st) .* _damped_sine_wave.(x, f, A, st, ϕ, d),
                            o .+ _ramp.(x, δ, st, et, h)], i)

@info "Testing voltage function generators..."
@testset "Voltage function generators" begin

    @named voltage_sensor = VoltageSensor()
    @named source = VoltageSource()

    @named step   = StepFunction(starttime=st, offset=o, height=h)
    @named square = SquareFunction(offset=o, starttime=st, amplitude=A, frequency=f)
    @named tri    = TriangularFunction(offset=o, starttime=st, amplitude=A, frequency=f)
    @named cosine = CosineFunction(offset=o, amplitude=A, frequency=f, starttime=st, phase=ϕ)
    @named ramp   = RampFunction(offset=o, starttime=st, endtime=et, height=h)
    @named damped_sine = DampedSineFunction(offset=o, amplitude=A, frequency=f, starttime=st, phase=ϕ, damping_coef=d)
    # @named vsawtooth = SawToothVoltage(amplitude=A, starttime=st, frequency=f, offset=o)
    
    wavefunctions = [tri, square, step, cosine, damped_sine, ramp]
    source_name = Symbolics.getname(source)
    @info "Testing for $source_name sources with following shapes of waves"
    for w in 1:length(wavefunctions)
        wave = wavefunctions[w]
        @info Symbolics.getname(wave)
        rc_eqs = [
            wave.y ~ source.v
            connect(voltage_sensor.p, source.p, resistor.p)
            connect(voltage_sensor.n, source.n, ground.g, resistor.n)
        ]
        @named model = ODESystem(rc_eqs, t, systems = [wave, source, resistor, voltage_sensor, ground])
        sys = structural_simplify(model)

        u0 = [
            source.v => 1
            resistor.v => 1
        ]

        prob = ODEProblem(sys, u0, (0, 10.0))
        sol = solve(prob, dt=0.1, Rosenbrock23())

        @test sol[source.v][1150:end] ≈ waveforms(w, sol.t)[1150:end] atol=1e-1
        # source_name == "current_source" && @test sol[source.i][1150:end] ≈ waveforms(w, sol.t)[1150:end] atol=1e-1
        # For visual inspection
        # plt = plot(sol)
        # savefig(plt, "test_current_$(Symbolics.getname(source))")
    end
end

@info "Testing current function generators..."
@testset "Current function generators" begin
    @named current_sensor = CurrentSensor()
    @named source = CurrentSource()

    @named step   = SmoothStepFunction(starttime=st, offset=o, height=h)
    @named square = SmoothSquareFunction(offset=o, starttime=st, amplitude=A, frequency=f)
    @named tri    = SmoothTriangularFunction(offset=o, starttime=st, amplitude=A, frequency=f)
    @named cosine = SmoothCosineFunction(offset=o, amplitude=A, frequency=f, starttime=st, phase=ϕ)
    @named ramp   = SmoothRampFunction(offset=o, starttime=st, endtime=et, height=h)
    @named damped_sine = SmoothDampedSineFunction(offset=o, amplitude=A, frequency=f, starttime=st, phase=ϕ, damping_coef=d)
    # @named vsawtooth = SawToothVoltage(amplitude=A, starttime=st, frequency=f, offset=o)
    
    wavefunctions = [tri, square, step, cosine, damped_sine, ramp]
    for w in 1:length(wavefunctions)
        wave = wavefunctions[w]
        @info "Testing $(Symbolics.getname(wave))"
        eqs = [
            wave.y ~ source.i
            connect(source.p, resistor.p)
            connect(source.n, resistor.n)
        ]
        @named model = ODESystem(eqs, t, systems = [wave, source, resistor])
        sys = alias_elimination(model)

        u0 = [
            source.i => 1
            resistor.v => 1
            source.v => 1
        ]

        prob = ODEProblem(sys, u0, (0, 10.0))
        sol = solve(prob, dt=0.1, Rosenbrock23())
        # plot(sol)
        @test sol[source.i][1150:end] ≈ waveforms(w, sol.t)[1150:end] atol=1e-1
    end
end