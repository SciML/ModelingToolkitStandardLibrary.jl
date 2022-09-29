using ModelingToolkit, OrdinaryDiffEq, Test

import ModelingToolkitStandardLibrary.Mechanical.Translational as TV
import ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition as TP

@parameters t
D = Differential(t)

@testset "spring damper mass fixed" begin

    @named dv = TV.Damper(d=1, v1₀=1)
    @named dp = TP.Damper(d=1, v1₀=1, s1₀=3, s2₀=1)

    @named sv = TV.Spring(k=1, v1₀=1, Δs₀=1)
    @named sp = TP.Spring(k=1, s1₀=3, s2₀=1, l=1)

    @named bv = TV.Mass(m=1, v₀=1)
    @named bp = TP.Mass(m=1, v₀=1, s₀=3)

    @named gv =  TV.Fixed()
    @named gp =  TP.Fixed(s₀=1)

    function simplify_and_solve(damping, spring, body, ground)

        eqs = [connect(spring.T1, body.T, damping.T1)
               connect(spring.T2, damping.T2, ground.T)
               ]
 
        @named model = ODESystem(eqs, t; systems = [ground, body, spring, damping])
 
        sys = structural_simplify(model)
 
        println.(full_equations(sys))
 
        prob = ODEProblem(sys, [], (0, 20.0), [])
        sol = solve(prob, ImplicitMidpoint(), dt=0.01)
 
        return sol
    end

    solv=simplify_and_solve(dv, sv, bv, gv);
    solp=simplify_and_solve(dp, sp, bp, gp);

    @test solv[bv.v][1] == 1.0
    @test solv[bv.v][end] ≈ 0.0 atol = 1e-4

    @test solp[bp.v][1] == 1.0
    @test solp[bp.v][end] ≈ 0.0 atol = 1e-4


end
