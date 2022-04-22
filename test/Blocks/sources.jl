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
    @test sol[src.output.u][end] ≈ 2 atol=1e-3
end

@testset "Sine" begin
    sine(t, frequency, amplitude, phase, offset, start_time) = offset + ifelse(t < start_time, 0, amplitude* sin(2*pi*frequency*(t - start_time) + phase))

    frequency=1
    amplitude=2
    phase=0
    offset=1
    start_time=0

    @named src = Sine(frequency=frequency, amplitude=amplitude, phase=phase, offset=offset, start_time=start_time)
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
    @test sol[src.output.u] ≈ sine.(sol.t, frequency, amplitude, phase, offset, start_time) atol=1e-3
end

@testset "Cosine" begin
    cosine(t, frequency, amplitude, phase, offset, start_time) = offset + ifelse(t < start_time, 0, amplitude* cos(2*pi*frequency*(t - start_time) + phase))
    
    frequency=1
    amplitude=2
    phase=0
    offset=1
    start_time=0

    @named src = Cosine(frequency=frequency, amplitude=amplitude, phase=phase, offset=offset, start_time=start_time)
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
    @test sol[src.output.u] ≈ cosine.(sol.t, frequency, amplitude, phase, offset, start_time) atol=1e-3
end

@testset "ContinuousClock" begin
    cont_clock(t, offset, start_time) = offset + ifelse(t < start_time, 0, t - start_time)
    
    offset, start_time = 1, 0

    @named src = ContinuousClock(offset=offset, start_time=start_time)
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
    @test sol[src.output.u] ≈ cont_clock.(sol.t, offset, start_time) atol=1e-3
end

@testset "Ramp" begin
    ramp(t, offset, height, duration, start_time) = offset + ifelse(t < start_time, 0, ifelse(t < (start_time + duration), (t - start_time) * height / duration, height))
    
    offset, height, duration, start_time = 1, 2, 2, 0

    @named src = Ramp(offset=offset, height=height, duration=duration, start_time=start_time)
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
    @test sol[src.output.u] ≈ ramp.(sol.t, offset, height, duration, start_time) atol=1e-3
end

@testset "Step" begin
    step(t, offset, height, start_time) = offset + ifelse(t < start_time, 0, height)

    offset, height, start_time = 1, 2, 5

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
    exp_sine(t, amplitude, frequency, damping, phase, start_time) = offset + ifelse(t < start_time, 0, amplitude * exp(-damping * (t - start_time)) * sin(2*pi*frequency*(t - start_time) + phase))

    frequency, amplitude, damping, phase, offset, start_time = 3, 2, 0.10, 0, 0, 0

    @named src = ExpSine(frequency=frequency, amplitude=amplitude, damping=damping, phase=phase, offset=offset, start_time=start_time)
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
    @test sol[src.output.u] ≈ exp_sine.(sol.t, amplitude, frequency, damping, phase, start_time) atol=1e-3
end