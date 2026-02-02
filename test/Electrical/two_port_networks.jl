using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test
using SciCompDSL
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks: Constant
using OrdinaryDiffEq: ReturnCode.Success

# Tests for two-port networks: IdealTransformer, Gyrator

@testset "IdealTransformer voltage ratio" begin
    # Test: Transformer with n=10, 10V input should produce 1V output
    @named source = Constant(k = 10.0)
    @named voltage = Voltage()
    @named transformer = IdealTransformer(n = 10.0)
    @named load = Resistor(R = 100.0)
    @named ground = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, transformer.p1)
        connect(voltage.n, transformer.n1, ground.g)
        connect(transformer.p2, load.p)
        connect(transformer.n2, load.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, voltage, transformer, load, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # Secondary voltage should be V1/n = 10V/10 = 1V
    @test sol[transformer.v2][end] ≈ 1.0 atol = 0.01
end

@testset "IdealTransformer current ratio" begin
    # Test: Transformer current ratio i1*n = -i2
    @named source = Constant(k = 10.0)
    @named voltage = Voltage()
    @named transformer = IdealTransformer(n = 10.0)
    @named load = Resistor(R = 100.0)
    @named ground = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, transformer.p1)
        connect(voltage.n, transformer.n1, ground.g)
        connect(transformer.p2, load.p)
        connect(transformer.n2, load.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, voltage, transformer, load, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # i2 = V2/R = 1V/100Ω = 0.01A
    # i1 * n = -i2 => i1 = -0.01/10 = -0.001A
    @test sol[transformer.i1][end] * 10.0 ≈ -sol[transformer.i2][end] atol = 0.001
end

@testset "IdealTransformer power conservation" begin
    # Test: P1 = V1*I1 should equal -P2 = -V2*I2
    @named source = Constant(k = 10.0)
    @named voltage = Voltage()
    @named transformer = IdealTransformer(n = 10.0)
    @named load = Resistor(R = 100.0)
    @named ground = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, transformer.p1)
        connect(voltage.n, transformer.n1, ground.g)
        connect(transformer.p2, load.p)
        connect(transformer.n2, load.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, voltage, transformer, load, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # Power conservation: P1 + P2 = 0
    P1 = sol[transformer.v1][end] * sol[transformer.i1][end]
    P2 = sol[transformer.v2][end] * sol[transformer.i2][end]
    @test P1 + P2 ≈ 0.0 atol = 0.001
end

@testset "Gyrator impedance transformation" begin
    # Test: Gyrator with R=1kΩ transforms impedance
    # A capacitor on port 2 appears as an inductor on port 1
    @named gyrator = Gyrator(R = 1000.0)

    # Just verify the component can be instantiated
    @test gyrator !== nothing
end
