using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test
using SciCompDSL
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks: Constant
using OrdinaryDiffEq: ReturnCode.Success

# Tests for controlled sources: VCVS, VCCS, CCVS, CCCS

@testset "VCVS voltage gain verification" begin
    # Test: VCVS with G=10, 1V input should produce 10V output
    @named source = Constant(k = 1.0)
    @named voltage = Voltage()
    @named vcvs = VCVS(G = 10.0)
    @named load = Resistor(R = 1000.0)
    @named ground = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, vcvs.p1)
        connect(voltage.n, vcvs.n1, ground.g)
        connect(vcvs.p2, load.p)
        connect(vcvs.n2, load.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, voltage, vcvs, load, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # Output voltage should be 10 * 1V = 10V
    @test sol[vcvs.v2][end] ≈ 10.0 atol = 0.01
    # Input current should be zero (infinite input impedance)
    @test sol[vcvs.i1][end] ≈ 0.0 atol = 1e-10
end

@testset "VCCS transconductance verification" begin
    # Test: VCCS with Gm=0.01 S, 2V input should produce 20mA output current
    @named source = Constant(k = 2.0)
    @named voltage = Voltage()
    @named vccs = VCCS(Gm = 0.01)
    @named load = Resistor(R = 100.0)
    @named ground = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, vccs.p1)
        connect(voltage.n, vccs.n1, ground.g)
        connect(vccs.p2, load.p)
        connect(vccs.n2, load.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, voltage, vccs, load, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # Output current should be 0.01 * 2V = 0.02A = 20mA
    @test sol[vccs.i2][end] ≈ 0.02 atol = 0.001
    # Input current should be zero
    @test sol[vccs.i1][end] ≈ 0.0 atol = 1e-10
end

@testset "CCVS transimpedance verification" begin
    # Test: CCVS with Rm=1000 Ω, 5mA input current should produce 5V output
    @named source = Constant(k = 0.005)  # 5mA current source
    @named current_src = Current()
    @named ccvs = CCVS(Rm = 1000.0)
    @named load = Resistor(R = 1000.0)
    @named ground = Ground()

    connections = [
        connect(source.output, current_src.I)
        connect(current_src.p, ccvs.p1)
        connect(current_src.n, ccvs.n1, ground.g)
        connect(ccvs.p2, load.p)
        connect(ccvs.n2, load.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, current_src, ccvs, load, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # Output voltage should be 1000 * 0.005A = 5V
    @test sol[ccvs.v2][end] ≈ 5.0 atol = 0.01
    # Input voltage should be zero (zero sensing impedance)
    @test sol[ccvs.v1][end] ≈ 0.0 atol = 1e-10
end

@testset "CCCS current gain verification" begin
    # Test: CCCS with α=100, 10μA input current should produce 1mA output
    @named source = Constant(k = 10e-6)  # 10μA current source
    @named current_src = Current()
    @named cccs = CCCS(α = 100.0)
    @named load = Resistor(R = 1000.0)
    @named ground = Ground()

    connections = [
        connect(source.output, current_src.I)
        connect(current_src.p, cccs.p1)
        connect(current_src.n, cccs.n1, ground.g)
        connect(cccs.p2, load.p)
        connect(cccs.n2, load.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, current_src, cccs, load, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # Output current should be 100 * 10μA = 1mA = 0.001A
    @test sol[cccs.i2][end] ≈ 0.001 atol = 1e-5
    # Input voltage should be zero (zero sensing impedance)
    @test sol[cccs.v1][end] ≈ 0.0 atol = 1e-10
end

@testset "Controlled sources with negative gain" begin
    # Test that controlled sources work with negative gain (inverting)
    @named source = Constant(k = 1.0)
    @named voltage = Voltage()
    @named vcvs = VCVS(G = -5.0)  # Inverting gain
    @named load = Resistor(R = 1000.0)
    @named ground = Ground()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, vcvs.p1)
        connect(voltage.n, vcvs.n1, ground.g)
        connect(vcvs.p2, load.p)
        connect(vcvs.n2, load.n, ground.g)
    ]

    @named model = System(connections, t;
        systems = [source, voltage, vcvs, load, ground])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, [], (0.0, 0.1))
    sol = solve(prob, Rodas4())

    @test SciMLBase.successful_retcode(sol)
    # Output voltage should be -5 * 1V = -5V
    @test sol[vcvs.v2][end] ≈ -5.0 atol = 0.01
end
