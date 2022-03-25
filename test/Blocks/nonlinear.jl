using ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks: Saturation, DeadZone, Integrator

@parameters t

#=
Testing strategy:
The general strategy is to test systems using simple intputs where the solution is known on closed form. For algebraic systems (without differential variables), an integrator with a constant input is often used together with the system under test. 
=#

@testset "Saturation" begin
    y_max = 0.8
    y_min = -0.6

    @named c = Constant(; k=1)
    @named int = Integrator(; k=1)
    @named sat = Saturation(; y_min, y_max)
    @named iosys = ODESystem([connect(c.y, int.u), connect(int.y, sat.u)], t, systems=[int, c, sat])
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, [int.x=>1.0], (0.0, 1.0))

    sol = solve(prob, Rodas4(), saveat=0:0.1:1)
    @test sol[int.y.u][end] ≈ 2
    @test sol[sat.y.u][end] ≈ 0.8
end

@testset "DeadZone" begin
    u_max = 1
    u_min = -2
    
    @named c = Constant(; k=1)
    @named int = Integrator(; k=1)
    @named dz = DeadZone(; u_min, u_max)
    @named iosys = ODESystem([connect(c.y, int.u), connect(int.y, dz.u)], t, systems=[int, c, dz])
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, [int.x=>1.0], (0.0, 1.0))

    sol = solve(prob, Rodas4(), saveat=0:0.1:1)
    @test sol[int.y.u][end] ≈ 2
end