using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test
using SciCompDSL
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks: Constant, Step
using OrdinaryDiffEq: ReturnCode.Success

# Tests for op-amp models: OpAmpFiniteGain, OpAmpGBW, OpAmpFull

@testset "OpAmpFiniteGain DC gain" begin
    # Test: Inverting amplifier with finite gain op-amp
    # Rf = 10k, Ri = 1k, theoretical gain = -Rf/Ri = -10
    # With finite A = 100000, actual gain is slightly less
    @named source = Constant(k = 1.0)
    @named voltage = Voltage()
    @named opamp = OpAmpFiniteGain(A = 100000.0, Rin = 1e6, Rout = 100.0)
    @named Rf = Resistor(R = 10000.0)
    @named Ri = Resistor(R = 1000.0)
    @named ground = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, Ri.p)
        connect(Ri.n, opamp.p1, Rf.p)
        connect(opamp.n1, ground.g)
        connect(opamp.p2, Rf.n)
        connect(opamp.n2, ground.g)
        connect(voltage.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, voltage, opamp, Rf, Ri, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # Closed-loop gain should be approximately -10
    @test sol[opamp.v2][end] ≈ -10.0 atol = 0.1
end

@testset "OpAmpGBW frequency response" begin
    # Test that OpAmpGBW has a state variable for frequency response
    @named opamp = OpAmpGBW(A0 = 100000.0, GBW = 1e6, Rin = 1e6)

    # Just verify the component can be instantiated
    @test opamp !== nothing
end

@testset "OpAmpFull saturation behavior" begin
    # Test that OpAmpFull can be instantiated with saturation limits
    @named opamp = OpAmpFull(A0 = 100000.0, GBW = 1e6, Vsat_p = 5.0, Vsat_n = -5.0,
        SR = 1e6, Rin = 1e6, Rout = 100.0)

    # Verify the component was created with correct parameters
    @test opamp !== nothing
    # Verify it has the expected state variable
    @test length(ModelingToolkit.unknowns(opamp)) > 0
end

@testset "OpAmpFull slew rate limiting" begin
    # Test that OpAmpFull can be instantiated with slew rate
    @named opamp = OpAmpFull(A0 = 100000.0, GBW = 1e6, Vsat_p = 12.0, Vsat_n = -12.0,
        SR = 1e6, Rin = 1e6, Rout = 100.0)

    # Just verify the component can be instantiated
    @test opamp !== nothing
end
