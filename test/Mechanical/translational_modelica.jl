using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

using ModelingToolkitStandardLibrary.Blocks: Sine
using ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica: Damper, Spring, Mass,
                                                                       Fixed, Force

@testset "spring damper mass fixed" begin
    @mtkmodel SpringDamperMassFixed begin
        @components begin
            damper = Damper(; d = 1)
            spring = Spring(; c = 1, s_rel0 = 1)
            mass = Mass(; m = 1, v = 1, s = 0)
            fixed = Fixed(s0 = 1)
        end
        @equations begin
            connect(spring.flange_a, mass.flange_a, damper.flange_a)
            connect(spring.flange_b, damper.flange_b, fixed.flange)
        end
    end

    @mtkbuild sys = SpringDamperMassFixed()

    prob = ODEProblem(sys, [], (0, 20.0), [])
    sol = solve(prob, ImplicitMidpoint(), dt = 0.01)

    @test sol[sys.mass.v][1] == 1.0
    @test sol[sys.mass.v][end]≈0.0 atol=1e-4
end

@testset "driven spring damper mass" begin
    @mtkmodel DrivenSpringDamperMass begin
        @components begin
            damper = Damper(; d = 1)
            spring = Spring(; c = 1, s_rel0 = 1)
            mass = Mass(; m = 1, v = 1, s = 0)
            fixed = Fixed(; s0 = 1)
            force = Force()
            source = Sine(frequency = 3, amplitude = 2)
        end

        @equations begin
            connect(force.f, source.output)
            connect(force.flange, mass.flange_a)
            connect(spring.flange_a, mass.flange_b, damper.flange_a)
            connect(spring.flange_b, damper.flange_b, fixed.flange)
        end
    end

    @mtkbuild sys = DrivenSpringDamperMass()

    prob = ODEProblem(sys, [], (0, 20.0), [])
    sol = solve(prob, Rodas4())

    lb, ub = extrema(sol(15:0.05:20, idxs = sys.mass.v).u)
    @test -lb≈ub atol=1e-2
    @test -0.11 < lb < -0.1
end
