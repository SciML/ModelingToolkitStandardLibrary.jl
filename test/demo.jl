using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq #, Plots
using Test

@testset "RC demo" begin
    R = 1.0
    C = 1.0
    V = 1.0
    @parameters t
    @named resistor = Resistor(R=R)
    @named capacitor = Capacitor(C=C)
    @named source = ConstantVoltage(V=V)
    @named ground = Ground()

    rc_eqs = [
            connect(source.p, resistor.p)
            connect(resistor.n, capacitor.p)
            connect(capacitor.n, source.n, ground.g)
            ]

    @named rc_model = ODESystem(rc_eqs, t, systems=[resistor, capacitor, source, ground])
    sys = structural_simplify(rc_model)
    prob = ODAEProblem(sys, Pair[], (0, 10.0))
    sol = solve(prob, Tsit5())
    #plot(sol)

    @test isapprox(sol[capacitor.v][end], V, atol=1e-2)
end