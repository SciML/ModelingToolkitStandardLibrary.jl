using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkitStandardLibrary.Blocks: _step, _square_wave, _triangular_wave,
                                 _cos_wave, _damped_sine_wave, _ramp
# using Plots

@parameters t
@named ground = Ground()
R = 5
@named resistor = Resistor(R=R)

@info "Testing the sensors..."
@testset "sensors" begin
    offset = 1
    freq = 50
    @named source = SineVoltage(offset=offset, amplitude=1.0, frequency=freq, starttime=.5, phase=0.0)

    @testset "current" begin
        @named current_sensor = CurrentSensor()
        
        rc_eqs = [
            connect(source.p, resistor.p)
            connect(resistor.n, current_sensor.p)
            connect(current_sensor.n, source.n, ground.g)
            ]
            @named rc_model = ODESystem(rc_eqs, t, systems = [resistor, source, current_sensor, ground])
            sys = structural_simplify(rc_model)
            u0 = [
                resistor.p.i => 10.0
                ]
                prob = ODEProblem(sys, u0, (0, 2.0))
                sol = solve(prob, Rosenbrock23())
                
                @test sol[current_sensor.i][1:20000:end] ≈ [(t < 0.5) ? offset/R : (offset + sin(2π*freq*(t .- 0.5)))/R for t in sol.t[1:20000:end]] atol=1e-6
            end
            
    @testset "potential" begin
        @named potential_sensor1 = PotentialSensor()
        @named potential_sensor2 = PotentialSensor()

        rc_eqs = [
                connect(source.p, resistor.p, potential_sensor1.p)
                connect(resistor.n, potential_sensor2.p, source.n, ground.g)
                ]
        @named rc_model = ODESystem(rc_eqs, t, systems = [resistor, source, potential_sensor1, potential_sensor2, ground])
        sys = structural_simplify(rc_model)
        u0 = [
            resistor.p.i => 10.0
        ]
        prob = ODEProblem(sys, u0, (0, 2.0))
        sol = solve(prob, Rosenbrock23())
        @test sol[potential_sensor1.phi][1:20000:end] ≈ [(t < 0.5) ? offset : offset + sin(2π*freq*(t .- 0.5)) for t in sol.t[1:20000:end]] atol=1e-6
        @test iszero(sol[potential_sensor2.phi][1:20000:end])
    end

    @testset "voltage" begin
        @named voltage_sensor = VoltageSensor()

        rc_eqs = [
                connect(source.p, resistor.p, voltage_sensor.p)
                connect(voltage_sensor.n, source.n, resistor.n, ground.g)
                ]
        @named rc_model = ODESystem(rc_eqs, t, systems = [resistor, source, voltage_sensor, ground])
        sys = structural_simplify(rc_model)
        u0 = [
            resistor.p.i => 10.0
        ]
        prob = ODEProblem(sys, u0, (0, 2.0))
        sol = solve(prob, Rosenbrock23())
        @test sol[voltage_sensor.v][1:20000:end] ≈ [(t < 0.5) ? offset : offset + sin(2π*freq*(t .- 0.5)) for t in sol.t[1:20000:end]] atol=1e-6
    end
    
    @testset "power" begin
        @named power_sensor = PowerSensor()

        rc_eqs = [
                connect(source.p, resistor.p, power_sensor.pv)
                connect(power_sensor.nv, resistor.n, power_sensor.nc)
                connect(power_sensor.pc, source.n, ground.g)
                ]
        @named rc_model = ODESystem(rc_eqs, t, systems = [resistor, source, power_sensor, ground])
        sys = structural_simplify(rc_model)
        u0 = [
            resistor.p.i => 10.0
        ]
        prob = ODEProblem(sys, u0, (0, 2.0))
        sol = solve(prob, Rosenbrock23())
        @test sol[power_sensor.power][1:20000:end] ≈ [(t < 0.5) ? offset^2/R : (offset + sin(2π*freq*(t .- 0.5)))^2 / R for t in sol.t[1:20000:end]] atol=1e-6
    end

    @testset "multi" begin
        @named multi_sensor = MultiSensor()

        rc_eqs = [
                connect(source.p, resistor.p, multi_sensor.pv)
                connect(multi_sensor.nv, resistor.n, multi_sensor.nc)
                connect(multi_sensor.pc, source.n, ground.g)
                ]
        @named rc_model = ODESystem(rc_eqs, t, systems = [resistor, source, multi_sensor, ground])
        sys = structural_simplify(rc_model)
        u0 = [
            resistor.p.i => 10.0
        ]
        prob = ODEProblem(sys, u0, (0, 2.0))
        sol = solve(prob, Rosenbrock23())
        @test sol[multi_sensor.i][1:20000:end] ≈ [(t < 0.5) ? offset/R : (offset + sin(2π*freq*(t .- 0.5)))/R for t in sol.t[1:20000:end]] atol=1e-6
        @test sol[multi_sensor.v][1:20000:end] ≈ [(t < 0.5) ? offset : offset + sin(2π*freq*(t .- 0.5)) for t in sol.t[1:20000:end]] atol=1e-6
    end
end

@info "Testing the inductor..."
@testset "Inductor" begin
    freq, offset = 5000, 1
    @named sinesource = SineVoltage(offset=offset, starttime=0.5, amplitude=1.0, frequency=5, phase=0.0)
    @named l1 = Inductor()
    @named l2 = Inductor()
    @named voltage_sensor = VoltageSensor()
    l_eqs = [
        connect(voltage_sensor.p, sinesource.p, l1.p)
        connect(l1.n, l2.p)
        connect(voltage_sensor.n, sinesource.n, l2.n, ground.g)
            ]
    @named l_model = ODESystem(l_eqs, t, systems = [l1, l2, sinesource, voltage_sensor, ground])
    sys = structural_simplify(l_model)
    equations(sys)
    u0 = [
        l1.p.i => 10.0
        sinesource.v => 1.0
        l2.v => 0.0
    ]
    prob = ODEProblem(sys, u0, (0, 10.0))
    sol = solve(prob, Rosenbrock23())
    
    @test sol[l1.v] + sol[l2.v] ≈ sol[voltage_sensor.v]
end

@info "Constructing an integrator circuit..."
# This tests Capacitor, IdealOpAmp
@testset "Integrator" begin
    @named ground = Ground()
    @named res1 = Resistor()
    @named c1 = Capacitor()
    @named opamp = IdealOpAmp()
    @named square = SquareVoltage()
    
    in_eqs = [
        connect(square.p, res1.p)
        connect(res1.n, c1.p, opamp.p1)
        connect(opamp.n2, c1.n)
        connect(opamp.n1, ground.g, opamp.p2, square.n)
    ]
    @named in_model = ODESystem(in_eqs, t, systems = [res1,  opamp, square, c1, ground])
    sys = structural_simplify(in_model)
    u0 = [
        c1.v => 1
        res1.v => 1
    ]
    prob = ODEProblem(sys, u0, (0, 10.0))
    sol = solve(prob, Rosenbrock23())
    @test sol[opamp.v2] == sol[c1.v] # Not a great one however. Rely on the plot

    # plt = plot(sol)
    # savefig(plt, "integrator")
end

#=@info "Testing the Current generators..."
@testset "Current function generators" begin
    st, o, h, f, A, et, ϕ, d, δ = 0.7, 1.25, 3, 2, 2.5, 2.5, π/4, 0.1, 0.0001

    @named ground = Ground()
    @named res = Resistor()
    @named current_sensor = CurrentSensor()
    @named istep = StepCurrent(starttime=st, offset=o, height=h)
    @named isquare = SquareCurrent(offset=o, starttime=st, amplitude=A, frequency=f)
    @named itri = TriangularCurrent(offset=o, starttime=st, amplitude=A, frequency=f)
    # @named isawtooth = SawToothCurrent(amplitude=A, starttime=st, frequency=f, offset=o)
    @named icosine = CosineCurrent(offset=o, amplitude=A, frequency=f, starttime=st, phase=ϕ)
    @named idamped_sine = DampedSineCurrent(offset=o, amplitude=A, frequency=f, starttime=st, phase=ϕ, damping_coef=d)
    @named iramp = RampCurrent(offset=o, starttime=st, endtime=et, height=h)

    isources = [itri, isquare, istep, icosine, idamped_sine, iramp]
    waveforms(i, x) = getindex([o .+ (x .> st) .* _triangular_wave.(x, δ, f, A, st),
                                o .+ (x .> st) .* _square_wave.(x, δ, f, A, st),
                                o .+  _step.(x, δ, h, st),
                                # o .+ (x .> st). * _sawtooth_wave.(x, δ, f, A, st),
                                o .+ (x .> st) .* _cos_wave.(x, f, A, st, ϕ),
                                o .+ (x .> st) .* _damped_sine_wave.(x, f, A, st, ϕ, d),
                                o .+ _ramp.(x, δ, st, et, h)], i)

    for i in 1:length(isources)
        isource = isources[i]
        eqs = [
            connect(isource.p, current_sensor.n)
            connect(current_sensor.p, res.n)
            connect(isource.n, res.p)
        ]
        @named model = ODESystem(eqs, t, systems = [current_sensor, isource, res])
        isys = alias_elimination(model)

        u0 = [
            isource.i => 1
            res.v => 1
        ]
        prob = ODEProblem(isys, u0, (0, 2.0))
        sol = solve(prob, Rosenbrock23())

        @test sol[isource.i][1150:end] ≈ waveforms(i, sol.t)[1150:end] atol=1e-1
        # For visual inspection
        # plt = plot(sol)
        # savefig(plt, "test_current_$(Symbolics.getname(source))")
    end
end=#