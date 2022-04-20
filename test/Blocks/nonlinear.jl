using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks

@parameters t

@testset "Limiter" begin
    y_max = 0.8
    y_min = -0.6

    @named c = Constant(; k=1)
    @named int = Integrator(; k=1)
    @named sat = Limiter(; y_min, y_max)
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

@testset "SlewRateLimiter" begin
    @named source = Sine(; frequency=1/2)
    @named rl = SlewRateLimiter(; rising=1, falling=-1, Td=0.001, y_start=-1/3)
    @named iosys = ODESystem([
        connect(source.output, rl.input),
        ],
        t,
        systems=[source, rl],
    )
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[], (0.0, 10.0))

    sol = solve(prob, Rodas4(), saveat=0.01, abstol=1e-10, reltol=1e-10)
    @test all(abs.(sol[rl.output.u]) .<= 0.51)
end