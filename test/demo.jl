using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq #, Plots
using ModelingToolkitStandardLibrary.Blocks: Constant
using Test

@testset "RC demo" begin
    R = 1.0
    C = 1.0
    V = 1.0
    @parameters t
    @named resistor  = Resistor(R=R)
    @named capacitor = Capacitor(C=C)
    @named voltage   = Voltage()
    @named ground    = Ground()
    @named source    = Constant()

    rc_eqs = [
            connect(source.output, voltage.V)
            connect(voltage.p, resistor.p)
            connect(resistor.n, capacitor.p)
            connect(capacitor.n, voltage.n, ground.g)
            ]

    @named rc_model = ODESystem(rc_eqs, t, systems=[resistor, capacitor, source, voltage, ground])
    sys = structural_simplify(rc_model)
    prob = ODAEProblem(sys, Pair[], (0, 10.0))
    sol = solve(prob, Tsit5())
    #plot(sol)
    @test sol.retcode == :Success
    @test isapprox(sol[capacitor.v][end], V, atol=1e-2)
end