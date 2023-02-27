using ModelingToolkit, OrdinaryDiffEq, Test

using ModelingToolkitStandardLibrary.Blocks
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as TP

@parameters t
D = Differential(t)

@testset "spring damper mass fixed" begin
    @named damper = TP.Damper(1)

    @named spring = TP.Spring(1; s_rel0 = 1)

    @named mass = TP.Mass(1)

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
    @named damper = TP.Damper(d = 1, v_a_0 = 1, s_a_0 = 3, s_b_0 = 1)

    @named spring = TP.Spring(k = 1, s_a_0 = 3, s_b_0 = 1, l = 1)

    @named mass = TP.Mass(m = 1, v_0 = 1, s_0 = 3)

    @named fixed = TP.Fixed(s_0 = 1)

    @named force = TP.Force()

    @named source = Sine(frequency = 3, amplitude = 2)

    function simplify_and_solve(damping, spring, body, fixed, f, source)
        eqs = [connect(f.f, source.output)
               connect(f.flange, body.flange)
               connect(spring.flange_a, body.flange, damping.flange_a)
               connect(spring.flange_b, damping.flange_b, fixed.flange)]

        @named model = ODESystem(eqs, t;
                                 systems = [fixed, body, spring, damping, f, source])

        sys = structural_simplify(model)

        foreach(println, full_equations(sys))

        prob = ODEProblem(sys, [], (0, 20.0), [])
        sol = solve(prob, Rodas4())

        return sol
    end

    solp = simplify_and_solve(damper, spring, mass, fixed, force, source)

    lb, ub = extrema(sol(15:0.05:20, idxs = mass.v).u)
    @test -lb≈ub atol=1e-2
    @test -0.11 < lb < -0.1

end
