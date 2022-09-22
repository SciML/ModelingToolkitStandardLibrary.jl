using ModelingToolkitStandardLibrary.Mechanical.Translational, ModelingToolkit,
      OrdinaryDiffEq,
      Test

@parameters t
D = Differential(t)

@testset "one falling body" begin
    #body starts at 1m, moving up at 2m/s in gravity of 10
    @named body = Body(x0 = 1.0, v0 = 2.0, g = -10.0)

    @named model = ODESystem(body.port.f ~ 0, t, [], []; systems = [body])

    sys = structural_simplify(model)
    prob = ODEProblem(sys, [], (0, 1.0))
    sol = solve(prob, ImplicitMidpoint(), dt = 0.1)

    #x = g*t^2/2 + v_int*t + x_int
    #dx = g*t + v_int

    @test sol.retcode == :Success
    @test sol[body.x, 1] == 1.0
    @test sol[body.dx, 1] == 2.0
    @test sol[body.x, end] ≈ -10 / 2 + 2.0 + 1.0
end
