using ModelingToolkit, OrdinaryDiffEq, Test

using ModelingToolkitStandardLibrary.Blocks
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as TP

@parameters t
D = Differential(t)

@testset "spring damper mass fixed" begin
    @named damper = TP.Damper(1)

    @named spring = TP.Spring(1; s_rel0 = 1)

    @named mass = TP.Mass(1, v0 = 1)

    @named fixed = TP.Fixed(s0 = 1)

    eqs = [connect(spring.flange_a, mass.flange_a, damper.flange_a)
           connect(spring.flange_b, damper.flange_b, fixed.flange)]

    @named model = ODESystem(eqs, t; systems = [fixed, mass, spring, damper])

    sys = structural_simplify(model)

    foreach(println, full_equations(sys))

    prob = ODEProblem(sys, [], (0, 20.0), [])
    sol = solve(prob, ImplicitMidpoint(), dt = 0.01)

    @test sol[mass.v][1] == 1.0
    @test sol[mass.v][end]≈0.0 atol=1e-4
end

@testset "driven spring damper mass" begin
    @named damper = TP.Damper(1)

    @named spring = TP.Spring(1; s_rel0 = 1)

    @named mass = TP.Mass(1, v0 = 1)

    @named fixed = TP.Fixed(s0 = 1)

    @named force = TP.Force()

    @named source = Sine(frequency = 3, amplitude = 2)

    eqs = [connect(force.f, source.output)
           connect(force.flange, mass.flange_a)
           connect(spring.flange_a, mass.flange_b, damper.flange_a)
           connect(spring.flange_b, damper.flange_b, fixed.flange)]

    @named model = ODESystem(eqs, t;
                             systems = [fixed, mass, spring, damper, force, source])

    sys = structural_simplify(model)

    foreach(println, full_equations(sys))

    prob = ODEProblem(sys, [], (0, 20.0), [])
    sol = solve(prob, Rodas4())

    lb, ub = extrema(sol(15:0.05:20, idxs = mass.v).u)
    @test -lb≈ub atol=1e-2
    @test -0.11 < lb < -0.1
end
