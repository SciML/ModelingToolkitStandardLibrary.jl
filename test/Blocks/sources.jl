using ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks

@parameters t

@testset "Constant" begin
    @named src = Constant(k=2)
    @named int = Integrator()
    @named iosys = ODESystem([
        connect(src.output, int.input),
        ],
        t,
        systems=[int, src],
    )
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol_lim[src.output.u][end] ≈ 2 atol=1e-3
end

@testset "Sine" begin
    @named src = Sine(frequency=1, amplitude=2, phase=0, offset=1, start_time=0)
    @named int = Integrator()
    @named iosys = ODESystem([
        connect(src.output, int.input),
        ],
        t,
        systems=[int, src],
    )
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol[src.output.u] ≈ 1 + 2 * sin(sol.t) atol=1e-3
end

@testset "Cosine" begin
    @named src = Cosine(frequency=1, amplitude=2, phase=0, offset=1, start_time=0)
    @named int = Integrator()
    @named iosys = ODESystem([
        connect(src.output, int.input),
        ],
        t,
        systems=[int, src],
    )
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol[src.output.u] ≈ 1 + 2 * cos(sol.t) atol=1e-3
end

@testset "ContinuousClock" begin
    @named src = ContinuousClock(offset=1, start_time=0)
    @named int = Integrator()
    @named iosys = ODESystem([
        connect(src.output, int.input),
        ],
        t,
        systems=[int, src],
    )
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol[src.output.u] ≈ 1 + sol.t atol=1e-3
end

@testset "Ramp" begin
    @named src = Ramp(offset=1, height=2, duration=2, start_time=0)
    @named int = Integrator()
    @named iosys = ODESystem([
        connect(src.output, int.input),
        ],
        t,
        systems=[int, src],
    )
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol[src.output.u][1] ≈ 1 atol=1e-3
    @test sol[src.output.u][end] ≈ 1 + 2 atol=1e-3
end

@testset "Step" begin
    step(t, offset, height, start_time) = offset + t < start_time ? 0 : height

    offset=1, height=2, start_time=5

    @named src = Step(offset=offset, height=height, start_time=start_time)
    @named int = Integrator()
    @named iosys = ODESystem([
        connect(src.output, int.input),
        ],
        t,
        systems=[int, src],
    )
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol[src.output.u] ≈ step.(sol.t, offset, height, start_time) atol=1e-3
end

@testset "ExpSine" begin
    @named src = ExpSine(frequency=3, amplitude=2, damping=0.1, phase=0, offset=0, start_time=0)
    @named int = Integrator()
    @named iosys = ODESystem([
        connect(src.output, int.input),
        ],
        t,
        systems=[int, src],
    )
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol[src.output.u] ≈ 2 * exp(-0.1*t) * sin(2*pi*3*t) atol=1e-3
end