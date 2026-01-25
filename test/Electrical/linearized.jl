using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test
using SciCompDSL
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks: Constant, Sine
using OrdinaryDiffEq: ReturnCode.Success

# Tests for linearized components: LinearizedDiode

@testset "LinearizedDiode forward bias behavior" begin
    # Test: Forward biased diode with Vd=0.7V, Rd=10Ω
    @named source = Constant(k = 2.0)  # 2V source
    @named voltage = Voltage()
    @named diode = LinearizedDiode(Vd = 0.7, Rd = 10.0)
    @named load = Resistor(R = 100.0)
    @named ground = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, diode.p)
        connect(diode.n, load.p)
        connect(load.n, voltage.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, voltage, diode, load, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # With 2V source, diode is forward biased
    # i = (V_source - Vd) / (Rd + Rload) = (2 - 0.7) / 110 ≈ 0.0118A
    @test sol[diode.i][end] > 0.0
end

@testset "LinearizedDiode reverse bias behavior" begin
    # Test: Reverse biased diode should have zero current
    @named source = Constant(k = -2.0)  # -2V source (reverse bias)
    @named voltage = Voltage()
    @named diode = LinearizedDiode(Vd = 0.7, Rd = 10.0)
    @named load = Resistor(R = 100.0)
    @named ground = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, diode.p)
        connect(diode.n, load.p)
        connect(load.n, voltage.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, voltage, diode, load, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # Reverse biased, current should be zero (or very small)
    @test sol[diode.i][end] ≈ 0.0 atol = 0.001
end

@testset "LinearizedDiode in half-wave rectifier" begin
    # Test that diode can be instantiated for rectifier circuits
    @named diode = LinearizedDiode(Vd = 0.7, Rd = 10.0)

    # Just verify the component can be instantiated
    @test diode !== nothing
end
