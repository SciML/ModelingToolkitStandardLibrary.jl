using ModelingToolkit, OrdinaryDiffEq, Test

using ModelingToolkitStandardLibrary.Blocks
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as TM
using DynamicQuantities: @u_str

@parameters t [unit = u"s"]
D = Differential(t)

@testset "spring damper mass fixed" begin
    @named damper = TM.Damper(; d = 1)
    @named spring = TM.Spring(; c = 1, s_rel0 = 1)
    @named mass = TM.Mass(; m = 1, v = 1)
    @named fixed = TM.Fixed(s0 = 1)

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
    @named damper = TM.Damper(; d = 1)
    @named spring = TM.Spring(; c = 1, s_rel0 = 1)
    @named mass = TM.Mass(; m = 1, v = 1)
    @named fixed = TM.Fixed(; s0 = 1)
    @named force = TM.Force(use_support = false)
    @named source = Sine(frequency = 3, amplitude = 2, output__unit = u"N")

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
