using ModelingToolkit, ModelingToolkitStandardLibrary, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Blocks: smooth_sin, smooth_cos, smooth_damped_sin,
                                             smooth_square, smooth_step, smooth_ramp,
                                             smooth_triangular, triangular, square
using OrdinaryDiffEq: ReturnCode.Success

@parameters t

@testset "Constant" begin
    @named src = Constant(k = 2)
    @named int = Integrator()
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u][end]≈2 atol=1e-3
end

@testset "TimeVaryingFunction" begin
    f(t) = t^2

    @named src = TimeVaryingFunction(f)
    @named int = Integrator()
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈f.(sol.t) atol=1e-3
    @test sol[int.output.u][end]≈1 / 3 * 10^3 atol=1e-3 # closed-form solution to integral
end

@testset "Sine" begin
    function sine(t, frequency, amplitude, phase, offset, start_time)
        offset + ifelse(t < start_time, 0,
               amplitude * sin(2 * pi * frequency * (t - start_time) + phase))
    end

    frequency = 1
    amplitude = 2
    phase = 0
    offset = 1
    start_time = 2
    δ = 1e-5
    @named int = Integrator()

    @named src = Sine(frequency = frequency, amplitude = amplitude, phase = phase,
                      offset = offset, start_time = start_time)
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈sine.(sol.t, frequency, amplitude, phase, offset, start_time) atol=1e-3

    @named smooth_src = Sine(frequency = frequency,
                             amplitude = amplitude,
                             phase = phase,
                             offset = offset,
                             start_time = start_time,
                             smooth = true)
    @named smooth_iosys = ODESystem([
                                        connect(smooth_src.output, int.input),
                                    ],
                                    t,
                                    systems = [int, smooth_src])

    smooth_sys = structural_simplify(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_sin.(smooth_sol.t, δ, frequency, amplitude,
                                                      phase, offset, start_time) atol=1e-3
end

@testset "Cosine" begin
    function cosine(t, frequency, amplitude, phase, offset, start_time)
        offset + ifelse(t < start_time, 0,
               amplitude * cos(2 * pi * frequency * (t - start_time) + phase))
    end

    frequency = 1
    amplitude = 2
    phase = 0
    offset = 1
    start_time = 2
    δ = 1e-5
    @named int = Integrator()

    @named src = Cosine(frequency = frequency,
                        amplitude = amplitude,
                        phase = phase,
                        offset = offset,
                        start_time = start_time,
                        smooth = false)
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])

    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈cosine.(sol.t, frequency, amplitude, phase, offset, start_time) atol=1e-3

    @named smooth_src = Cosine(frequency = frequency,
                               amplitude = amplitude,
                               phase = phase,
                               offset = offset,
                               start_time = start_time,
                               smooth = true)
    @named smooth_iosys = ODESystem([
                                        connect(smooth_src.output, int.input),
                                    ],
                                    t,
                                    systems = [int, smooth_src])

    smooth_sys = structural_simplify(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_cos.(smooth_sol.t, δ, frequency, amplitude,
                                                      phase, offset, start_time) atol=1e-3
end

@testset "ContinuousClock" begin
    cont_clock(t, offset, start_time) = offset + ifelse(t < start_time, 0, t - start_time)

    offset, start_time = 1, 0

    @named src = ContinuousClock(offset = offset, start_time = start_time)
    @named int = Integrator()
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])
    sys = structural_simplify(iosys)

    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))

    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈cont_clock.(sol.t, offset, start_time) atol=1e-3
end

@testset "Ramp" begin
    function ramp(t, offset, height, duration, start_time)
        offset + ifelse(t < start_time, 0,
               ifelse(t < (start_time + duration), (t - start_time) * height / duration,
                      height))
    end

    offset, height, duration, start_time, δ = 1, 2, 2, 0, 1e-5
    @named int = Integrator()

    @named src = Ramp(offset = offset, height = height, duration = duration,
                      start_time = start_time)
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])

    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈ramp.(sol.t, offset, height, duration, start_time) atol=1e-3

    start_time = 2
    @named smooth_src = Ramp(offset = offset, height = height, duration = duration,
                             start_time = start_time, smooth = true)
    @named smooth_iosys = ODESystem([
                                        connect(smooth_src.output, int.input),
                                    ],
                                    t,
                                    systems = [int, smooth_src])

    smooth_sys = structural_simplify(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_ramp.(smooth_sol.t, δ, height, duration,
                                                       offset, start_time) atol=1e-3
end

@testset "Step" begin
    step(t, offset, height, start_time) = offset + ifelse(t < start_time, 0, height)

    offset, height, start_time, δ = 1, 2, 5, 1e-5
    @named int = Integrator()

    @named src = Step(offset = offset, height = height, start_time = start_time,
                      smooth = false)
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])

    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())

    @test sol.retcode == Success
    @test sol[src.output.u]≈step.(sol.t, offset, height, start_time) atol=1e-2

    # test with duration
    duration = 1.2
    @named src = Step(offset = offset, height = height, start_time = start_time,
                      duration = duration, smooth = false)
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])

    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4(), dtmax = 0.1) # set dtmax to prevent the solver from overstepping the entire step disturbance

    @test sol.retcode == Success
    @test sol[src.output.u]≈step.(sol.t, offset, height, start_time) -
                            step.(sol.t, 0, height, start_time + duration) atol=1e-2

    @named smooth_src = Step(offset = offset, height = height, start_time = start_time,
                             smooth = true)
    @named smooth_iosys = ODESystem([
                                        connect(smooth_src.output, int.input),
                                    ],
                                    t,
                                    systems = [int, smooth_src])

    smooth_sys = structural_simplify(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4(), dtmax = 0.1) # set dtmax to prevent the solver from overstepping the entire step disturbance)

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u] ≈
          smooth_step.(smooth_sol.t, δ, height, offset, start_time)

    # with duration
    @named smooth_src = Step(offset = offset, height = height, start_time = start_time,
                             smooth = true, duration = duration)
    @named smooth_iosys = ODESystem([
                                        connect(smooth_src.output, int.input),
                                    ],
                                    t,
                                    systems = [int, smooth_src])

    smooth_sys = structural_simplify(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u] ≈
          smooth_step.(smooth_sol.t, δ, height, offset, start_time) -
          smooth_step.(smooth_sol.t, δ, height, 0, start_time + duration)
end

@testset "Square" begin
    frequency = 1
    amplitude = 2
    offset = 1
    start_time = 2.5
    δ = 1e-5
    @named int = Integrator()

    @named src = Square(frequency = frequency, amplitude = amplitude,
                        offset = offset, start_time = start_time)
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])

    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())

    @test sol.retcode == Success
    @test sol[src.output.u]≈square.(sol.t, frequency, amplitude, offset, start_time) atol=1e-3

    @named smooth_src = Square(frequency = frequency, amplitude = amplitude,
                               offset = offset, start_time = start_time, smooth = true)
    @named smooth_iosys = ODESystem([
                                        connect(smooth_src.output, int.input),
                                    ],
                                    t,
                                    systems = [int, smooth_src])

    smooth_sys = structural_simplify(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_square.(smooth_sol.t, δ, frequency,
                                                         amplitude, offset, start_time) atol=1e-3
end

@testset "Triangular" begin
    frequency = 5
    amplitude = 1
    offset = 2
    start_time = 1
    δ = 1e-5
    @named int = Integrator()

    @named src = Triangular(frequency = frequency, amplitude = amplitude,
                            offset = offset, start_time = start_time)
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])

    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 4.0))
    sol = solve(prob, Rodas4(), saveat = 0.01)

    @test sol.retcode == Success
    @test sol[src.output.u]≈triangular.(sol.t, frequency, amplitude, offset, start_time) atol=1e-3

    @named smooth_src = Triangular(frequency = frequency, amplitude = amplitude,
                                   offset = offset, start_time = start_time, smooth = true)
    @named smooth_iosys = ODESystem([
                                        connect(smooth_src.output, int.input),
                                    ],
                                    t,
                                    systems = [int, smooth_src])

    smooth_sys = structural_simplify(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 4.0))
    smooth_sol = solve(smooth_prob, Rodas4(), saveat = 0.01)

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_triangular.(smooth_sol.t, δ, frequency,
                                                             amplitude, offset, start_time) atol=1e-3
end

@testset "ExpSine" begin
    function exp_sine(t, amplitude, frequency, damping, phase, start_time)
        offset + ifelse(t < start_time, 0,
               amplitude * exp(-damping * (t - start_time)) *
               sin(2 * pi * frequency * (t - start_time) + phase))
    end

    frequency, amplitude, damping = 3, 2, 0.10
    phase, offset, start_time, δ = 0, 0, 0, 1e-5
    @named src = ExpSine(frequency = frequency, amplitude = amplitude, damping = damping,
                         phase = phase, offset = offset, start_time = start_time)
    @named int = Integrator()
    @named iosys = ODESystem([
                                 connect(src.output, int.input),
                             ],
                             t,
                             systems = [int, src])
    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[int.x => 0.0], (0.0, 10.0))
    sol = solve(prob, Rodas4())
    @test sol.retcode == Success
    @test sol[src.output.u]≈exp_sine.(sol.t, amplitude, frequency, damping, phase,
                                      start_time) atol=1e-3

    offset, start_time = 1, 2
    @named smooth_src = ExpSine(frequency = frequency, amplitude = amplitude,
                                damping = damping, phase = phase, offset = offset,
                                start_time = start_time, smooth = true)
    @named smooth_iosys = ODESystem([
                                        connect(smooth_src.output, int.input),
                                    ],
                                    t,
                                    systems = [int, smooth_src])
    smooth_sys = structural_simplify(smooth_iosys)
    smooth_prob = ODEProblem(smooth_sys, Pair[int.x => 0.0], (0.0, 10.0))
    smooth_sol = solve(smooth_prob, Rodas4())

    @test smooth_sol.retcode == Success
    @test smooth_sol[smooth_src.output.u]≈smooth_damped_sin.(smooth_sol.t, δ, frequency,
                                                             amplitude, damping, phase,
                                                             offset, start_time) atol=1e-3
end
