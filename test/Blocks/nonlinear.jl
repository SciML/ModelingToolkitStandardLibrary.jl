using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks

@parameters t

@testset "Saturation" begin
    y_max = 0.8
    y_min = -0.6

    @named c = Constant(; k=1)
    @named int = Integrator(; k=1)
    @named sat = Saturation(; y_min, y_max)
    @named model = ODESystem([
            connect(c.output, int.input), 
            connect(int.output, sat.input),
        ], 
        t, 
        systems=[int, c, sat],
    )
    sys = structural_simplify(model)

    prob = ODEProblem(sys, [int.x=>1.0], (0.0, 1.0))

    sol = solve(prob, Rodas4())
    @test sol[int.output.u][end] ≈ 2
    @test sol[sat.output.u][end] ≈ 0.8
end

@testset "DeadZone" begin
    u_max = 1
    u_min = -2
    
    @named c = Constant(; k=1)
    @named int = Integrator(; k=1)
    @named dz = DeadZone(; u_min, u_max)
    @named model = ODESystem([
            connect(c.output, int.input), 
            connect(int.output, dz.input),
        ], 
        t, 
        systems=[int, c, dz],
    )
    sys = structural_simplify(model)

    prob = ODEProblem(sys, [int.x=>1.0], (0.0, 1.0))

    sol = solve(prob, Rodas4())
    @test sol[int.output.u][end] ≈ 2
end