using ModelingToolkit, OrdinaryDiffEq, Test

using ModelingToolkitStandardLibrary.Blocks
import ModelingToolkitStandardLibrary: Mechanical
import ModelingToolkitStandardLibrary.Mechanical.Translational as TV
import ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition as TP
using DynamicQuantities: @u_str

@parameters t [unit = u"s"]
D = Differential(t)

@testset "Free" begin
    function System(; name)
        systems = @named begin
            acc = TV.Acceleration()
            a = Constant(; k = -10, output.unit = u"m/s^2")
            mass = TV.Mass(; m = 100)
            free = TV.Free()
        end

        eqs = [connect(a.output, acc.a)
            connect(mass.flange, acc.flange, free.flange)]

        ODESystem(eqs, t, [], []; name, systems)
    end

    @named system = System()
    s = complete(system)
    sys = structural_simplify(system)
    prob = ODEProblem(sys, [], (0, 0.1))
    sol = solve(prob, Rosenbrock23())

    @test sol[s.mass.flange.v][end]≈-0.1 * 10 atol=1e-3
    @test sol[s.free.f][end] ≈ 100 * 10
end

@testset "Spring, Damper, Mass, Fixed" begin
    @named dv = TV.Damper(d = 1, flange_a.v = 1)
    @named dp = TP.Damper(d = 1, va = 1, vb = 0.0, flange_a.s = 3, flange_b.s = 1)

    @named sv = TV.Spring(k = 1, flange_a__v = 1, delta_s = 1)
    @named sp = TP.Spring(k = 1, flange_a__s = 3, flange_b__s = 1, l = 1)

    @named bv = TV.Mass(m = 1, v = 1)
    @named bp = TP.Mass(m = 1, v = 1, s = 3)

    @named gv = TV.Fixed()
    @named gp = TP.Fixed(s_0 = 1)

    function simplify_and_solve(damping, spring, body, ground)
        eqs = [connect(spring.flange_a, body.flange, damping.flange_a)
            connect(spring.flange_b, damping.flange_b, ground.flange)]

        @named model = ODESystem(eqs, t; systems = [ground, body, spring, damping])

        sys = structural_simplify(model)

        prob = ODEProblem(sys, [], (0, 20.0), [])
        sol = solve(prob, ImplicitMidpoint(), dt = 0.01)

        return sol
    end

    solv = simplify_and_solve(dv, sv, bv, gv)
    solp = simplify_and_solve(dp, sp, bp, gp)

    @test solv[bv.v][1] == 1.0
    @test solv[bv.v][end]≈0.0 atol=1e-4

    @test solp[bp.v][1] == 1.0
    @test solp[bp.v][end]≈0.0 atol=1e-4
end

@testset "driven spring damper mass" begin
    @named dv = TV.Damper(d = 1, flange_a.v = 1)
    @named dp = TP.Damper(d = 1, va = 1.0, vb = 0.0, flange_a.s = 3, flange_b.s = 1)

    @named sv = TV.Spring(k = 1, flange_a__v = 1, delta_s = 1)
    @named sp = TP.Spring(k = 1, flange_a__s = 3, flange_b__s = 1, l = 1)

    @named bv = TV.Mass(m = 1, v = 1)
    @named bp = TP.Mass(m = 1, v = 1, s = 3)

    @named gv = TV.Fixed()
    @named gp = TP.Fixed(s_0 = 1)

    @named fv = TV.Force()
    @named fp = TP.Force(use_support = false)

    @named source = Sine(frequency = 3, amplitude = 2, unit = u"N")

    function System(damping, spring, body, ground, f, source)
        eqs = [connect(f.f, source.output)
            connect(f.flange, body.flange)
            connect(spring.flange_a, body.flange, damping.flange_a)
            connect(spring.flange_b, damping.flange_b, ground.flange)]

        @named model = ODESystem(eqs, t;
            systems = [ground, body, spring, damping, f, source])

        return model
    end

    model = System(dv, sv, bv, gv, fv, source)
    sys = structural_simplify(model)
    prob = ODEProblem(sys, [], (0, 20.0), [])
    solv = solve(prob, Rodas4())

    model = System(dp, sp, bp, gp, fp, source)
    sys = structural_simplify(model)
    prob = ODEProblem(sys, [], (0, 20.0), [])
    solp = solve(prob, Rodas4())

    for sol in (solv, solp)
        lb, ub = extrema(solv(15:0.05:20, idxs = bv.v).u)
        @test -lb≈ub atol=1e-2
        @test -0.11 < lb < -0.1
    end
end

@testset "sources & sensors" begin
    function System(; name)
        systems = @named begin
            pos = TV.Position()
            pos_sensor = TV.PositionSensor(; s = 1)
            force = TV.Force()
            force_sensor = TV.ForceSensor()

            spring = TV.Spring(; k = 1000)

            src1 = Sine(frequency = 100, amplitude = 2, output__unit = u"m")
            src2 = Sine(frequency = 100, amplitude = -1, output__unit = u"N")

            pos_value = RealInput(unit = u"m")
            force_output = RealOutput(unit = u"N")
        end

        eqs = [connect(pos.s, src1.output)
            connect(force.f, src2.output)
            connect(spring.flange_a, pos.flange, force_sensor.flange)
            connect(spring.flange_b, force.flange, pos_sensor.flange)
            connect(pos_value, pos_sensor.output)
            connect(force_output, force_sensor.output)]

        ODESystem(eqs, t, [], []; name, systems)
    end

    @named system = System()
    s = complete(system)
    sys = structural_simplify(system)
    prob = ODEProblem(sys, [], (0, 1 / 400))
    sol = solve(prob, Rosenbrock23())

    delta_s = 1 / 1000
    s_b = 2 - delta_s + 1

    @test sol[s.pos_value.u][end]≈s_b atol=1e-3
end
