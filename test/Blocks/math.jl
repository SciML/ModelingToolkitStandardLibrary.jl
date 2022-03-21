using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks

@parameters t

#=
Testing strategy:
The general strategy is to test systems using simple intputs where the solution is known on closed form. For algebraic systems (without differential variables), an integrator with a constant input is often used together with the system under test. 
=#

@testset "Gain" begin
    @info "Testing Gain"
    @named c = Gain(42)
    @named int = Integrator(; k=1)
    @named iosys = ODESystem([int.u~c.y, c.u~1], t, systems=[c, int])
    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[], (0.0, 1.0))
    sol = solve(prob, Rosenbrock23(), saveat=0:0.1:1)
    @test sol[int.y] ≈ 42 .* (0:0.1:1)

    # Matrix gain
    @named c = Gain([2 0; 0 3])
    ints = [Integrator(; k=1, name=Symbol("int$i")) for i in 1:2]
    @named iosys = ODESystem([
        ints[1].u~c.y[1],
        ints[2].u~c.y[2],
        c.u[1]~1,
        c.u[2]~2,
        ], t, systems=[c; ints])
    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[], (0.0, 1.0))
    sol = solve(prob, Rosenbrock23(), saveat=0:0.1:1)
    @test sol[ints[1].y] ≈ 2 .* (0:0.1:1) # 2 * 1
    @test sol[ints[2].y] ≈ 6 .* (0:0.1:1) # 3 * 2
end


@testset "Sum" begin
    @info "Testing Sum"
    @named s = Sum(2)
    ints = [Integrator(; k=1, name=Symbol("int$i")) for i in 1:2]
    @named iosys = ODESystem([
        ints[1].u~1,
        ints[2].u~2,
        ints[1].y~s.u[1],
        ints[2].y~s.u[2],
        ], t, systems=[s; ints])
    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[], (0.0, 1.0))
    sol = solve(prob, Rosenbrock23(), saveat=0:0.1:1)
    @test sol[s.y] ≈ 3 .* (0:0.1:1)

    @named s = Sum([1, -2])
    ints = [Integrator(; k=1, name=Symbol("int$i")) for i in 1:2]
    @named iosys = ODESystem([
        ints[1].u~1,
        ints[2].u~1,
        ints[1].y~s.u[1],
        ints[2].y~s.u[2],
        ], t, systems=[s; ints])
    sys = structural_simplify(iosys)
    prob = ODEProblem(sys, Pair[], (0.0, 1.0))
    sol = solve(prob, Rosenbrock23(), saveat=0:0.1:1)
    @test sol[s.y] ≈ (1 + (-2)) .* (0:0.1:1)
end