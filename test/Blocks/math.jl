using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit, OrdinaryDiffEq

@parameters t


@testset "Gain" begin
    @named c = Constant(; k=1)
    @named gain = Gain(1;)
    @named int = Integrator(; k=1)
    @named model = ODESystem([connect(c.output, gain.input), connect(gain.output, int.input)], t, systems=[int, gain, c])
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[int.x=>1.0], (0.0, 1.0))

    sol = solve(prob, Rodas4())

    @test sol[int.output.u][end] ≈ 2
end

@testset "Feedback loop" begin
    @named c = Constant(; k=2)
    @named gain = Gain(1;)
    @named int = Integrator(; k=1)
    @named fb = Feedback(;)
    @named model = ODESystem(
        [
            connect(c.output, fb.input1), 
            connect(fb.input2, int.output), 
            connect(fb.output, gain.input),
            connect(gain.output, int.input),
        ], 
        t, 
        systems=[int, gain, c, fb]
    )
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 100.0))

    sol = solve(prob, Rodas4())
    @test sol[int.output.u][end] ≈ 2
end

@testset "Add" begin
    @named c1 = Constant(; k=1)
    @named c2 = Constant(; k=2)
    @named add = Add(;)
    @named int = Integrator(; k=1)
    @named model = ODESystem(
        [
            connect(c1.output, add.input1), 
            connect(c2.output, add.input2), 
            connect(add.output, int.input),
        ], 
        t, 
        systems=[int, add, c1, c2]
    )
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 1.0))

    sol = solve(prob, Rodas4())
    @test sol[int.output.u][end] ≈ 3
end

@testset "Product" begin
    @named c1 = Constant(; k=1)
    @named c2 = Constant(; k=2)
    @named prod = Product(;)
    @named int = Integrator(; k=1)
    @named model = ODESystem(
        [
            connect(c1.output, prod.input1), 
            connect(c2.output, prod.input2), 
            connect(prod.output, int.input),
        ], 
        t, 
        systems=[int, prod, c1, c2]
    )
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 1.0))

    sol = solve(prob, Rodas4())
    @test sol[int.output.u][end] ≈ 2
end

@testset "Division" begin
    @named c1 = Constant(; k=1)
    @named c2 = Constant(; k=2)
    @named div = Division(;)
    @named int = Integrator(; k=1)
    @named model = ODESystem(
        [
            connect(c1.output, div.input1), 
            connect(c2.output, div.input2), 
            connect(div.output, int.input),
        ], 
        t, 
        systems=[int, div, c1, c2]
    )
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 1.0))

    sol = solve(prob, Rodas4())
    @test sol[int.output.u][end] ≈ 1/2
end

@testset "Abs" begin
    @named c = Constant(; k=-1)
    @named abs = Abs(;)
    @named int = Integrator(; k=1)
    @named model = ODESystem(
        [
            connect(c.output, abs.input), 
            connect(abs.output, int.input),
        ], 
        t, 
        systems=[int, abs, c]
    )
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 1.0))

    sol = solve(prob, Rodas4())
    @test sol[int.output.u][end] ≈ 1
end

@testset "Sqrt" begin
    @named c = Constant(; k=4)
    @named sqr = Sqrt(;)
    @named int = Integrator(; k=1)
    @named model = ODESystem(
        [
            connect(c.output, sqr.input), 
            connect(sqr.output, int.input),
        ], 
        t, 
        systems=[int, sqr, c]
    )
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 1.0))

    sol = solve(prob, Rodas4())
    @test sol[int.output.u][end] ≈ 2
end

@testset "Sign" begin
    @named c = Constant(; k=3)
    @named sig = Sign(;)
    @named int = Integrator(; k=1)
    @named model = ODESystem(
        [
            connect(c.output, sig.input), 
            connect(sig.output, int.input),
        ], 
        t, 
        systems=[int, sig, c]
    )
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 1.0))

    sol = solve(prob, Rodas4())
    @test sol[int.output.u][end] ≈ 1
end

@testset "MatrixGain" begin
    K = [1 2; 3 4]
    @named gain = MatrixGain(K;)
    # TODO:
end

@testset "Sum" begin
    @named s = Sum(2;)
    # TODO:
end

@testset "Math" begin
    blocks = [Abs, Sign, Sqrt, Sin, Cos, Tan, Asin, Acos, Atan, Sinh, Cosh, Tanh, Exp]
    for block in blocks
        @named source = Sine()
        @named b = block()
        @named int = Integrator()
        @named model = ODESystem([connect(source.output, b.input), connect(b.output, int.input)], t, systems=[int, b, source])
        sys = structural_simplify(model)

        prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 1.0))

        @test_nowarn sol = solve(prob, Rodas4())
    end

    blocks = [Log, Log10] # input must be positive
    for block in blocks
        @named source = Sin(; offset=2)
        @named b = block()
        @named int = Integrator()
        @named model = ODESystem([connect(source.output, b.input), connect(b.output, int.input)], t, systems=[int, b, source])
        sys = structural_simplify(model)

        prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 1.0))

        @test_nowarn sol = solve(prob, Rodas4())
    end
end

@testset "Atan2" begin
    @named c1 = Constant(; k=1)
    @named c2 = Constant(; k=2)
    @named b = Atan2(;)
    @named int = Integrator(; k=1)
    @named model = ODESystem(
        [
            connect(c1.output, b.input1), 
            connect(c2.output, b.input2), 
            connect(b.output, int.input),
        ], 
        t, 
        systems=[int, b, c1, c2]
    )
    sys = structural_simplify(model)

    prob = ODEProblem(sys, Pair[int.x=>0.0], (0.0, 1.0))

    sol = solve(prob, Rodas4())
    @test sol[int.output.u][end] ≈ atan(1, 2)
end